class FeatureError(Exception):
    def __init__(self,message="Feature not implemented"):
        self.message = message
        super().__init__(self.message)

class MyDockingError(Exception):
    def __init__(self,message="Error has not been assigned"):
        self.message = "!!! Error !!!!\n\t"+message
        super().__init__(self.message)
