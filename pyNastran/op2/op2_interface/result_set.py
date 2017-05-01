from copy import deepcopy

class ResultSet(object):
    """
    This class is private storage interface class.

    It's an interface tool between the code and the results the user requests.
    """
    def __init__(self, allowed_results):
        # the full set of allowable results
        self.allowed = set(allowed_results)

        # the set of results that have been found
        self.found = set([])

        # the set of results to be saved
        self.saved = deepcopy(self.allowed)

    def is_saved(self, result):
        if result not in self.allowed:
            #print(self.allowed)
            raise RuntimeError("result=%r is invalid; the name changed or it's a typo" % result)
        if result in self.saved:
            return True
        return False

    def is_not_saved(self, result):
        return not self.is_saved(result)

    def clear(self):
        self.saved.clear()

    def add(self, result):
        self.saved.add(result)

    def remove(self, results):
        for result in results:
            if result not in self.allowed:
                #allowed = list(self.allowed)
                #allowed.sort()
                #raise RuntimeError('%r is not a valid result to remove\nallowed=[%s]' % (
                    #result, ', '.join(allowed)))
                raise RuntimeError('%r is not a valid result to remove\n%s' % (
                                   result, self))

        for result in results:
            if result in self.saved:
                self.saved.remove(result)
        #disable_set = set(results)
        #self.saved.difference(disable_set)

    def _found_result(self, result):
        if result not in self.allowed:
            raise RuntimeError("result=%r is invalid; the name changed or it's a typo" % result)
        self.found.add(result)

    def update(self, results):
        for result in results:
            self.add(result)

    #def add_found_result(self, result):
        #pass

    def __repr__(self):
        msg = 'ResultSet:\n'
        msg += ' results:\n'
        for result in sorted(self.allowed):
            if result in self.saved:
                msg += '  %s\n' % result
            else:
                msg += '  %s (disabled)\n' % result
        return msg
