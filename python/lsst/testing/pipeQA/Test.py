
class Test(object):

    def __init__(self, label, value, limits, comment):
        self.label = label
        self.value = value
        self.limits = limits
        self.comment = comment


    def evaluate(self):
        """Add a test to this testing suite."""
        
        # grab a traceback for failed tests
        if (self.value < self.limits[0] or self.value > self.limits[1]):
            return False
        else:
            return True

