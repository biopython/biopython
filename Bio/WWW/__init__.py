# This is a Python module

class RequestLimiter:
    # This class implements a simple countdown timer for delaying WWW
    # requests.
    def __init__(self, delay):
        self.last_time = 0.0
        self.delay = delay
    def wait(self):
        how_long = self.last_time + self.delay - time.time()
        if how_long > 0:
            time.sleep(how_long)
        self.last_time = time.time()

