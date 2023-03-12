from nwalign import NWAligner
from swalign import SWAligner
import random

z = NWAligner('tgccacttaaaaaa', 'aaagtcagccc')
z.getAlignment()

x = SWAligner('tgccacttaaaaaa', 'aaagtcagccc')
x.getAlignment()