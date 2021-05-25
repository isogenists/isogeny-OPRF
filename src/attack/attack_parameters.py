from sage.all_cmdline import Integer

class AttackParameters():

    def __init__(self, eM = 4, eK = 3, e1 = 2, e2 = 3, e3 = 2, ec = 1, **kwargs):
        self.eM = Integer(eM)
        self.eK = Integer(eK)
        self.e1 = Integer(e1)
        self.e2 = Integer(e2)
        self.e3 = Integer(e3)
        self.ec = Integer(ec)
 
