import pylab
import modeller

def get_profile(profile_file, seq):
    """Read `profile_file` into a Python array, and add gaps corresponding to the alignment sequence `seq`."""
    #CG# read all non-comment and non-blank lines from the file:
    f = file(profile_file)
    vals = []
    for line in f:
        if not line.startswith('#') and len(line) > 10:
            spl = line.split()
            vals.append(float(spl[-1]))
    #CG# insert gaps into the profile corresponding to those in seq:
    for n, res in enumerate(seq.residues):
        for gap in range(res.get_leading_gaps()):
            vals.insert(n, None)
    #CG# add a gap at position '0', so that we effectively count from 1:
    vals.insert(0, None)
    return vals

e = modeller.environ()
a = modeller.alignment(e, file='res_align.ali')

template = get_profile('4GRV.profile', a['4GRV'])
model = get_profile('orexin.profile', a['orexin'])

#CG# plot the template and model profiles in the same plot for comparison:
pylab.figure(1, figsize=(10,6))
pylab.xlabel('Alignment position')
pylab.ylabel('DOPE per-residue score')
pylab.plot(model, color='red', linewidth=2, label='Model')
pylab.plot(template, color='green', linewidth=2, label='Template')
pylab.legend()
pylab.savefig('dope_profile.png', dpi=65)