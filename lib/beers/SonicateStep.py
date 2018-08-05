class SonicateStep:

    def __init__(self, log_file, parameters):
        print("Sonicate step instantiated")

    def execute(self, sample):
        print("Sonicate step acting on sample")
        return f"Post Sonicate step sample"

    def validate(self):
        print(f"Sonicate step validating parameters")
        return True
