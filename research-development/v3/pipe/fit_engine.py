import numpy as np
from numpy.linalg import svd, eig
from typing import Dict, Any

class AlignmentFitEngine:
    """
    Automatic global fitting engine for the Alignment Spectral Triple v2.0.

    Fits:
        • Quark masses
        • Lepton masses
        • CKM mixing matrix
        • PMNS mixing matrix
        • Neutrino hierarchy (optional)
    """

    def __init__(self, geom, triadic_kernel_cls, finite_cls, compressor):
        """
        Parameters
        ----------
        geom : GeometricTriple
        triadic_kernel_cls : class
            Reference to TriadicKernel class.
        finite_cls : class
            Reference to FiniteTriple class.
        compressor : TriadicCompression
        """
        self.geom = geom
        self.TriadicKernel = triadic_kernel_cls
        self.FiniteTriple = finite_cls
        self.compressor = compressor

    # ------------------------------------------------------------
    # Build all Yukawa matrices for a given parameter set
    # ------------------------------------------------------------
    def build_yukawas(self, params: Dict[str, Any]) -> Dict[str, np.ndarray]:
        """
        params includes:
            'sites'  : list of 9 ints in Z_2160
            'kappa'  : float
            'phases' : length-9 phase vector
            'R_u','R_d','R_e','R_nu' : 9-element weight vectors
        """
        # 1. Build real kernel K(sites, kappa)
        TK = self.TriadicKernel(
            kappa=params['kappa'],
            sites=params['sites'],
            forbidden=(2, 4, 7),
            N=2160
        )
        K_real = TK.K

        # 2. Phase-dressed D_F
        FT = self.FiniteTriple(K_real, params['phases'])
        D_F = FT.D_F

        # 3. Sector kernels (R K R)
        Ku  = self.compressor.sector_kernel(D_F, params['R_u'])
        Kd  = self.compressor.sector_kernel(D_F, params['R_d'])
        Ke  = self.compressor.sector_kernel(D_F, params['R_e'])
        Kn  = self.compressor.sector_kernel(D_F, params['R_nu'])

        # 4. Compress to 3×3 Yukawas
        Yu  = self.compressor.compress(Ku)
        Yd  = self.compressor.compress(Kd)
        Ye  = self.compressor.compress(Ke)
        Yn  = self.compressor.compress(Kn)

        return {'Yu': Yu, 'Yd': Yd, 'Ye': Ye, 'Yn': Yn}

    # ------------------------------------------------------------
    # Diagonalize Yukawas (SVD)
    # ------------------------------------------------------------
    @staticmethod
    def diagonalize(Y):
        """
        Y = U diag(m_i) V†
        masses = singular values (sorted)
        """
        U, s, Vh = svd(Y)
        return U, s, Vh.conjugate().T

    # ------------------------------------------------------------
    # Build CKM and PMNS
    # ------------------------------------------------------------
    def mixing_matrices(self, Ys: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        Uu, mu, _ = self.diagonalize(Ys['Yu'])
        Ud, md, _ = self.diagonalize(Ys['Yd'])

        Ue, me, _ = self.diagonalize(Ys['Ye'])
        Un, mn, _ = self.diagonalize(Ys['Yn'])

        Vckm = Uu.conjugate().T @ Ud
        Upmns = Ue.conjugate().T @ Un

        return dict(
            Vckm=Vckm,
            Upmns=Upmns,
            mu=mu,
            md=md,
            me=me,
            mn=mn
        )

    def extract_angles(self, U):
        """
        Numerically stable extraction of mixing angles.
        Works for CKM and PMNS.
        """
        s13 = np.clip(abs(U[0, 2]), 0, 1)
        c13 = np.sqrt(max(1 - s13 ** 2, 1e-15))

        s12 = np.clip(abs(U[0, 1]) / max(c13, 1e-15), 0, 1)
        s23 = np.clip(abs(U[1, 2]) / max(c13, 1e-15), 0, 1)

        return np.arcsin(s12), np.arcsin(s23), np.arcsin(s13)

    # ------------------------------------------------------------
    # χ² cost function
    # ------------------------------------------------------------
    def chi2(self, Ys, target):

        # 1. Yukawa -> mixing matrices
        try:
            result = self.mixing_matrices(Ys)
        except Exception:
            return 1e12

        # 2. Guard numerical stability
        for key in ['mu', 'md', 'me', 'mn']:
            if np.any(~np.isfinite(result[key])) or np.any(result[key] <= 0):
                return 1e12

        if np.any(~np.isfinite(result['Vckm'])):
            return 1e12
        if np.any(~np.isfinite(result['Upmns'])):
            return 1e12

        # 3. Extract angles robustly
        θ12, θ23, θ13 = self.extract_angles(result['Vckm'])
        t12, t23, t13 = self.extract_angles(result['Upmns'])

        # 4. Check stability
        if not np.isfinite(θ12 + θ23 + θ13 + t12 + t23 + t13):
            return 1e12

        # 5. Mass-part chi²
        chi = 0.0

        def add(pred, tgt):
            return (np.log(pred / tgt)) ** 2

        mu1, mu2, mu3 = result['mu']
        md1, md2, md3 = result['md']
        me1, mmu, mtau = result['me']

        chi += add(mu1, target['mu'])
        chi += add(mu2, target['mc'])
        chi += add(mu3, target['mt'])

        chi += add(md1, target['md'])
        chi += add(md2, target['ms'])
        chi += add(md3, target['mb'])

        chi += add(me1, target['me'])
        chi += add(mmu, target['mmu'])
        chi += add(mtau, target['mtau'])

        # 6. Angle chi²
        chi += (θ12 - target['ckm_12']) ** 2
        chi += (θ23 - target['ckm_23']) ** 2
        chi += (θ13 - target['ckm_13']) ** 2

        chi += (t12 - target['pmns_12']) ** 2
        chi += (t23 - target['pmns_23']) ** 2
        chi += (t13 - target['pmns_13']) ** 2

        return float(chi)

