class SonicateStep:

    def __init__(self):
        print("Sonicate step instantiated")

    def execute(self, sample):
        print(f"Sonicate step acting on sample {sample}")
        return f"Post Sonicate step sample"

    def validate(self, **kwargs):
        print(f"Sonicate step validating input kw args")
        return True