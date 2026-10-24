# R2-D10-literature — inputs for the lab's D10 gate (modular fixed points)

Item #22, literature-11. Its board notes (seq 238, plus the typo fix seq 239) are saved verbatim below by the lead,
2026-10-23. All reads are OPENED-img (PDF page images) unless marked. Item [5] was not opened by the member; the
lead opened it, see the addendum.

## Note (seq 238)

#22 D10 literature. Method: WebFetch saved the PDFs, I Read them as page images (OPENED-img). curl/pdftotext are refused for this role.
CORRECTION: arXiv:1907.11714 is Ding-King-Liu 'Modular A4 Symmetry Models of Neutrinos and Charged Leptons' (title p1), plain modular A4, NO gCP, no fixed-point analysis. Not usable for the ask. Fixed-point paper = Ding-King-Liu-Lu arXiv:1910.03460 'Modular S4 and A4 Symmetries and Their Fixed Points' (also no CP: Sec 6 'without CP').
[1] NPPT 1905.11970 (modular S4 + gCP), OPENED-img pp1-3,12-18. (a) p12 eq4.19: CP-conserving tau = boundary of D and imaginary axis; 'at tau_C=i the symmetry is enhanced to Z2^CP x Z2^S ... at the cusps Z2^CP x Z3^ST / Z2^CP x Z3^TS ... at tau_T=i inf Z2^CP x Z_N^T'. p12 eq4.17-4.22: at those tau 'resulting mass matrices are real ... CPV phases trivial, sin delta = sin a21 = sin a31 = 0' (0 or pi). (b) p14: 'only complex parameter ... g'/g'; gCP => Im(g'/g)=0 eq5.4; p18: 7 real parameters (vd alpha, beta/alpha, gamma/alpha, vu^2 g^2/Lambda, g'/g, Im tau, Re tau-phase), 12 observables. (c) The paper does NOT evaluate at a fixed point; best fit is tau = +-0.09922+1.016i (near i, Table 1 p15). sin^2th23 = 0.4862 (2sig 0.4848-0.4873; 3sig 0.484-0.4882) i.e. th23 about 44.2 deg NOT 45, tied to tau by correlations (Fig 2 p17, 'accurate predictions at a few percent'). (d) delta/pi = +-1.641 (3sig 1.627-1.656), a21/pi=+-0.3464, a31/pi=+-1.254, NO, m1=0.01211 eV, sum m=0.07833 eV. (e) fit data set/NuFIT version: NOT seen on pages I read (not stated p15-18); unknown. Fit N sigma=1.012.
[2] Ding-King-Liu-Lu 1910.03460 (S4, A4; no gCP; THREE different moduli tau_e, tau_atm, tau_sol each at a fixed point), OPENED-img pp1-3,16-18,22-24. (a) Table 3 p17 (S4): case A Ge=Z3^ST, tau_atm=i, tau_sol=2+omega; B Ge=Z3^ST, tau_atm=2+i, tau_sol=i; C Ge=Z4^(ST^2), tau_atm=-1+i, tau_sol=(9+i sqrt3)/14. p16: 'no phenomenologically viable models with one common fixed point tau_f in both neutrino and charged lepton sectors' (common value => four or six zeros, p23). (b) 3 real params ma, ms, eta (p24) with m1=0. (c) th23 NOT fixed: Table 4 p18 sin^2th23 = 0.591 (A), 0.553 (B), 0.596 (C); depends on eta, r=ms/ma (p18 'contour plots ... in plane r vs eta'). A,B give TM1: first column (sqrt(2/3), -1/sqrt6, -1/sqrt6) eq80. (d) Table 4: delta_CP/pi = -0.369 (A), -0.426 (B), -0.365 (C); beta/pi -0.564,-0.473,-0.608; delta is not +-90. (e) Table 4: sin^2th13 0.0227/0.0223/0.0226, sin^2th12 0.318/0.318/0.346; chi2min 3.804/0.74?: read 3.804, 0.609, 11.188 (A,B,C columns as printed: 3.804, 0.609, 11.188). NuFIT version not seen on pages read. A4 weight-6 example eq101: sin^2th23=0.563, delta=-0.399pi, with Y alignment not fixed uniquely (p22).
[3] Feruglio 1706.08749 (A4, no gCP), OPENED-img pp1,17-20,30. Not at a fixed point: Ex.2 tau=-0.195+1.0636i (eq45) sin^2th23=0.331, delta/pi=1.47, NO (eq46); Ex.3 (no flavons, one modulus) tau=0.008+0.98i (near i) sin^2th23=0.749, sin^2th13=0.001, delta/pi=1.25, IO (eq53-55), p20 'excluded by many standard deviations'. p20 sec3.2: tau=i has SO(2) invariance of the SL(2,R) action. Data set not seen.
[4] Stabilisation NPP 2201.02020, OPENED-img pp1-9,17-18. Potential: N=1 sugra, K=-L_K log(2 Im tau) eq2.4, W=L_W^3 H(tau)/eta^(2n) eq2.6, H=(j-1728)^(m/2) j^(n/3) P(j) eq2.7, V eq2.9/2.12 (n=3). Parameters: integers m,n and polynomial P (taken =1), scale L_V; p8: 'we take the simplest choice P(j)=1'. Results p17-18: (m,n)=(0,0),(1,1),(0,3) global minima at tau=1.2i, +-0.24+0.97i, i (CP-conserving, Ref [69]); n=0,m>0 give degenerate CP-breaking minima tau_min and -conj(tau_min) near left cusp omega: |u|min = 0.0145/(m+0.0625) with u=(tau-omega)/(tau-omega^2), arg(u_min) about -2pi/9 = -40 deg, 'independently of m'. m=2: |u|min=0.00705, epsilon about 0.0197, inside 1sig of the S4' model of their Ref [20]. So minima are NEAR omega (not at it), CP-breaking. p7 fn11: 'we do not identify their minimum values with the cosmological constant' => no dS/AdS statement on pages I read; V=e^K(...-3|W|^2) allows either sign, not discussed there. p17: results independent of the finite group/weights/reps. Intro p3-4 (descriptions of others, via this paper): Ref [58] flux-landscape T^6/(Z2xZ2') minima CP-conserving, cluster at omega; Ref [63] sugra with dilaton gives CP-violating minima near i (Ref numbers only, I did not open those).
[5] Kobayashi-Otsuka: NOT opened. Search hit only: 'Landscape of modular symmetric flavor models' Ishiguro-Kobayashi-Otsuka, JHEP Mar 2021 (doaj snippet, SUMM, arXiv id not confirmed by me; I believe 2011.09154 but did not open it). Not quoted. Open if you want it.
Gaps: no model found with gCP evaluated EXACTLY at a fixed point; at those points gCP forces trivial phases (NPPT p12). NuFIT versions of all three fits unread.

## Typo fix (seq 239)

#22 typo fix in [2](e): Table 4 p18 chi2min = 3.804 (A), 0.609?? NO: read columns as printed: chi2min A 3.804, B 0.609 is the ms/ma? Re-read: header order is eta/pi, ma, ms/ma, chi2min; row A: 0.818, 2.166, 1.722, 3.804; B: 0.758, 31.355, 0.074, 0.609; C: 0.533, 4.352, 7.334, 11.188. So chi2min = 3.804/0.609/11.188 (A/B/C); ignore the stray '0.74?'. Also the Table 4 sin^2th13/th12/th23 for A,B,C: 0.0227/0.318/0.591, 0.0223/0.318/0.553, 0.0226/0.346/0.596.

## Lead addendum: [5] Kobayashi–Otsuka (OPENED-img by the lead, 2026-10-23)

K. Ishiguro, T. Kobayashi, H. Otsuka, "Landscape of Modular Symmetric Flavor Models", arXiv:2011.09154v1
(18 Nov 2020). Title and authors were read from the PDF p.1. The journal ref (JHEP, 2021) came only via the
member's search snippet (SUMM) and is not confirmed.

- Setting (pdftotext, OPENED-text): Type IIB flux compactification on T^6/(Z2 x Z2'). Three-form flux superpotential
  (their Eq. 20) with the Kähler potential Eq. 17; 4D N=1 supergravity scalar potential. The parameters are integer
  flux quanta, with the tadpole N_flux bounded by N_flux^max. The modulus tau and the axio-dilaton S are stabilised
  together. This is a statistical landscape count of stable vacua, not a single potential with continuous parameters.
- Abstract p.1 (OPENED-text): "distributions of complex structure modulus tau determining the flavor structure are
  clustered at a fixed point with the residual Z3 symmetry in the SL(2,Z) fundamental region. Also, they are
  clustered at other specific points such as intersecting points between |tau|^2 = k/2 and Re tau = 0, +-1/4, +-1/2,
  although their probabilities are less than the Z3 fixed point. In general, CP-breaking vacua in the complex
  structure modulus are statistically disfavored in the string landscape. Among CP-breaking vacua, the values
  Re tau = +-1/4 are most favorable in particular when the axio-dilaton S is stabilized at Re S = +-1/4."
- Table 2, p.14 (OPENED-img), N_flux^max = 10, probability of stable vacua: tau = omega (-1/2, sqrt3/2) 62.3%;
  (0, sqrt3) 7.55; (-1/2, sqrt15/2) 7.55; (-1/4, sqrt15/4) 7.55; tau = i (0, 1) 5.66; then 1.89 each for
  (-1/2, 3sqrt3/2), (0, sqrt(3/2)), (0, sqrt6), (0, sqrt2), (-1/2, sqrt7/2).
- Table 3, p.14 (OPENED-img), N_flux^max = 1000: tau = omega 40.3%; (0, sqrt3) 7.55; (-1/2, sqrt15/2) 4.85;
  (-1/4, sqrt15/4) 4.85; tau = i 3.79; (-1/2, sqrt7/2) 2.43; (0, sqrt6), (0, sqrt(3/2)), (0, sqrt2) 1.88 each;
  (-1/2, 3sqrt3/2) 1.49.
- p.14 text: "the stable vacua are clustered at the Z3 fixed point, when N_flux^max is small ... other favorable
  modulus values such as the Z2 fixed point (Re tau, Im tau) = (0, 1) and (0, k) with k = 2, 3, 6. All of these are
  on boundary of the GL(2,Z) fundamental region."
- dS/AdS: not read.

## Lead's reading for the gate (no new claim; just the quotes combined)

1. No model was found with gCP evaluated *exactly* at a fixed point. At gCP-conserving tau, NPPT p.12 gives real mass
   matrices and trivial phases (sin delta = 0). So a gCP model sitting exactly at tau = i or omega predicts
   delta in {0, 180}, not +-90.
2. Where theta23 was quoted, it is not fixed by the residual symmetry:
   - NPPT best fit near i: sin^2 theta23 = 0.4862, tied to tau through correlations.
   - DKLL fixed-point models: 0.553–0.596, depending on eta and r.
   - Feruglio: near-i examples, 0.331 and 0.749.
3. Stabilisation:
   - NPP 2022: CP-breaking minima *near* omega, at |u| = 0.0145/(m + 0.0625), with CP-conserving minima at i for
     (m,n) = (0,3).
   - Ishiguro–Kobayashi–Otsuka: the flux landscape clusters at omega (40–62%), with tau = i at 3.8–5.7%.
4. Open: the NuFIT versions used by the three fits were not seen on the pages read.
