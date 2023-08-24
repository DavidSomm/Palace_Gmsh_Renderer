
class PALACE_Model:
    
    #create internal class variable
    _engine = None

    def __init__(self, model_name):
        self.model_name = model_name

    @staticmethod
    def init_engine():
        pass    