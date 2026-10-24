#!/usr/bin/env python3
"""tex_static_check.py: no LaTeX toolchain on this machine. Expands \\input from
theta23_obstruction_paper.tex and checks that every \\ref/\\eqref has exactly one \\label, every \\cite
has exactly one \\bibitem, and the braces balance. Run from the repo root."""
import re,collections
def strip(s): return re.sub(r'(?<!\\)%.*','',s)
def expand(f,seen=()):
    s=strip(open(f).read())
    return re.sub(r'\\input\{([^}]+)\}',lambda m: expand(m.group(1)+('' if m.group(1).endswith('.tex') else '.tex')),s)
t=expand('theta23_obstruction_paper.tex')
labels=collections.Counter(re.findall(r'\\label\{([^}]+)\}',t))
refs=set(k for c in re.findall(r'\\(?:eq|c|auto)?ref\{([^}]+)\}',t) for k in c.split(','))
bib=collections.Counter(re.findall(r'\\bibitem(?:\[[^]]*\])?\{([^}]+)\}',t))
cites=set(k.strip() for c in re.findall(r'\\cite[pt]?(?:\[[^]]*\])?\{([^}]+)\}',t) for k in c.split(','))
print('inputs expanded from theta23_obstruction_paper.tex:',re.findall(r'\\input\{([^}]+)\}',strip(open('theta23_obstruction_paper.tex').read())))
print(f'labels {len(labels)}, refs {len(refs)}, bibitems {len(bib)}, cite keys {len(cites)}')
print('refs with no label:',sorted(refs-set(labels)) or 'none')
print('labels defined more than once:',sorted(k for k,v in labels.items() if v>1) or 'none')
print('cites with no bibitem:',sorted(cites-set(bib)) or 'none')
print('bibitems defined more than once:',sorted(k for k,v in bib.items() if v>1) or 'none')
print('bibitems never cited:',sorted(set(bib)-cites) or 'none')
b=t.replace('\\{','').replace('\\}',''); print('brace balance (expanded):',b.count('{')-b.count('}'))
