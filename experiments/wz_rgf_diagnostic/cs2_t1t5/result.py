class Results:
    def __init__(self, stage):
        self.stage, self.items = stage, []

    def add(self, name, ok, cls=None, detail=""):
        self.items.append(dict(name=name, ok=bool(ok), cls=None if ok else cls, detail=str(detail)))
        return ok

    @property
    def failures(self):
        return [x for x in self.items if not x["ok"]]

    def summary(self):
        return dict(stage=self.stage, total=len(self.items), failed=len(self.failures))
