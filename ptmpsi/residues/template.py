import numpy as np
from ptmpsi.exceptions import MyDockingError

class Template:
    def __init__(self):
        self.name = None
        self.coordinates = None
        self.elements = None
        self.fullname = None
        self.charges = None
        self.backbone = np.zeros(3)
        self.nattach = np.zeros(3)
        self.cattach = np.zeros(3)
        self.pka = None
        self.chi1 = None
        self.chi2 = None

    def find(self,atom):
        if self.elements is None:
            raise MyDockingError("Template was not initialized")
        pos = next((idx for idx,val in np.ndenumerate(self.elements[:,1]) if val==atom),None)
        if pos is None:
            raise MyDockingError("Atom '{}' was not found in Template '{}'".format(atom,self.name))
        return pos[0]


    def find_coord(self,atom):
        return self.coordinates[self.find(atom)]

