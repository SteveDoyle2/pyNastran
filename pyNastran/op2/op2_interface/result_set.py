"""
Defines:
 - ResultSet(allowed_results)

   Attributes
   ----------
    - allowed
    - found
    - saved

   Methods
   -------
   - is_saved(result)
   - is_not_saved(result)
   - clear()
   - add(result)
   - remove(results)
   - _found_result(result)
   - update(self, results)

"""
import re
from copy import deepcopy


class ResultSet:
    """
    This class is private storage interface class.

    It's an interface tool between the code and the results the user requests.

    """
    def __init__(self, allowed_results, results_map, unused_log):
        #self.log = log
        #allowed_results.sort()
        #for a in allowed_results:
            #print(a)

        # the full set of allowable results
        self.allowed = set(allowed_results)
        #assert 'responses.convergence_data' in allowed_results

        # the set of results that have been found
        self.found = set()

        # the set of results to be saved
        self.saved = deepcopy(self.allowed)
        self.results_map = results_map

    def is_saved(self, result):
        """checks to see if a result is saved"""
        if result not in self.allowed:
            #allowed2 = list(self.allowed)
            #allowed2.sort()
            msg = "result=%r is invalid; the name changed or it's a typo.\n" % result
            if '.' in result:
                base, end = result.split('.', 1)
                #print(base, end)
                #print(self.allowed)
                #print(f'base={base} end={end}')
                #print(self.results_map)
                if base in self.results_map:
                    results_obj = self.results_map[base]
                    msg += 'Potential results include:\n - ' + '\n - '.join(results_obj.get_table_types())
                    assert result in results_obj.get_table_types()
                    #print(results_obj.get_table_types())
                raise RuntimeError(msg.rstrip())
        if result in self.saved:
            #self.log.debug('    %s is being read' % result)
            return True
        #self.log.debug('    %s was skipped' % result)
        return False

    def is_not_saved(self, result):
        """checks to see if a result is saved"""
        return not self.is_saved(result)

    def clear(self):
        """clears all the results"""
        self.saved.clear()

    def add(self, results):
        """addds a list/str of results"""
        all_matched_results = self._get_matched_results(results)
        for result in all_matched_results:
            if result not in self.saved:
                self.saved.add(result)

    def remove(self, results):
        """removes a list/str of results"""
        all_matched_results = self._get_matched_results(results)
        for result in all_matched_results:
            if result in self.saved:
                self.saved.remove(result)
        #disable_set = set(results)
        #self.saved.difference(disable_set)

    def _get_matched_results(self, results):
        """handles expansion of regexs"""
        if isinstance(results, str):
            results = [results]
        all_matched_results = []
        for result in results:
            # tack on a word boundary if we have a * at the beginning of the regex
            resulti = r'\w' + result if result.startswith('*') else result
            regex = re.compile(resulti)
            matched_results = list(filter(regex.match, self.allowed))
            if len(matched_results) == 0:
                #allowed = list(self.allowed)
                #allowed.sort()
                #raise RuntimeError('%r is not a valid result to remove\nallowed=[%s]' % (
                    #result, ', '.join(allowed)))
                raise RuntimeError(f'{result!r} is not a valid result to remove\n{self}\n'
                                   f'{result!r} is not a valid result to remove')
            all_matched_results.extend(matched_results)
        return all_matched_results

    def _found_result(self, result):
        if result not in self.allowed:
            msg = "result=%r is invalid; the name changed or it's a typo" % result
            raise RuntimeError(msg) # check line ~640 in op2_f06_common.py if this is a new result
        self.found.add(result)

    def update(self, results):
        for result in results:
            self.add(result)

    #def add_found_result(self, result):
        #pass

    def __repr__(self):
        """defines the repr"""
        msg = 'ResultSet:\n'
        msg += ' results:\n'
        for result in sorted(self.allowed):
            if result in self.saved:
                msg += '  %s\n' % result
            else:
                msg += '  %s (disabled)\n' % result
        return msg
