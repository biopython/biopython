"""Deal with various biological databases and services on the web.
"""
import time

class RequestLimiter:
    # This class implements a simple countdown timer for delaying WWW
    # requests.
    def __init__(self, delay):
        self.last_time = 0.0
        self.delay = delay
    def wait(self, delay=None):
        if delay is None:
            delay = self.delay
        how_long = self.last_time + delay - time.time()
        if how_long > 0:
            time.sleep(how_long)
        self.last_time = time.time()
