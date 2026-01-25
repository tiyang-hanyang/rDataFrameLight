class GeneralJob:
    registry = {} 

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        GeneralJob.registry[cls.__name__] = cls

    def __init__(self):
        self.periods = []
        self.datasets = []
        self.workflow = []
        self.outDir = []

    def declare(self):
        raise NotImplementedError