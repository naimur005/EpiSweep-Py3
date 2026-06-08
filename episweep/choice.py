# This has to be cleaned up

class Choice :
    # superclass for what goes at a position (AA / tube)
    pass

class SeqChoice (Choice) :
    # a sequence-level choice, with possible structural-level sub-choices (default: self)
    def __init__(self):
        self.choices = [self]
        
    def is_wt(self):
        return False

class AAChoice (SeqChoice) :
    def __init__(self, aa, wt):
        self.aa = aa; self.wt = wt
        self.aas = [aa]

    def is_wt(self):
        return self.wt

    def __str__(self):
        return self.aa

class TubeChoice (SeqChoice) :
    def __init__(self, degen, aas, wt):
        self.aa = degen; self.idx = 0
        self.aas = aas; self.wt = wt
        self.name = degen

    def is_wt(self):
        return self.wt

    def __str__(self):
        return self.aa

class AAChoice_idx (SeqChoice) :
    def __init__(self, aa, idx):
        self.aa = aa; self.idx = idx
        self.aas = [aa]
        self.name = aa + str(idx)

    def __str__(self):
        return self.name
