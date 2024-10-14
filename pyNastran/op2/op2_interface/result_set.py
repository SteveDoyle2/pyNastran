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
from typing import Any
from cpylog import SimpleLogger


class ResultSet:
    """
    This class is private storage interface class.

    It's an interface tool between the code and the results the user requests.

    """
    def __init__(self, allowed_results: list[str],
                 results_map: dict[str, Any],
                 log: SimpleLogger):
        """we'll set the allowed results

        Parameters
        ----------
        allowed_results : list[str]
            all results that are allowed
        results_map : dict[str, Any]
            combined objects (e.g., ato, psd, stress, force, ...)
        log : SimpleLogger
            the logging object

        """
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
        #self.log = log

    def is_saved(self, result: str) -> bool:
        """checks to see if a result is saved"""
        #assert result in self.results_map, f'result={result}'
        #print(self.saved)
        #self.saved = list(self.saved)
        #self.saved.sort()
        if result not in self.allowed:
            #allowed2 = list(self.allowed)
            #allowed2.sort()
            msg = f"result={result!r} is invalid; the name changed or it's a typo.\n"
            if '.' in result:
                base, unused_end = result.split('.', 1)
                #print(base, end)
                #print(self.allowed)
                #print(f'base={base} end={end}')
                #print(self.results_map)
                if base in self.results_map:
                    results_obj = self.results_map[base]
                    msg += 'Potential results include:\n - ' + '\n - '.join(results_obj.get_table_types())
                    assert result in results_obj.get_table_types(), result
                    #print(results_obj.get_table_types())
                raise RuntimeError(msg.rstrip())
        if result in self.saved:
            #log.debug('    %s is being read' % result)
            #self.log.warning('    %s is being read' % result)
            return True
        #log.debug('    %s was skipped' % result)
        #self.log.warning('    %s was skipped' % result)
        return False

    def is_not_saved(self, result: str) -> bool:
        """checks to see if a result is saved"""
        return not self.is_saved(result)

    def clear(self) -> None:
        """clears all the results"""
        #self.log.warning('clearing')
        self.saved.clear()
        #self.log.warning(f'saved={self.saved}')

    def add(self, results: str | list[str])  -> list[str]:
        """adds a list/str of results"""
        added = []
        if len(results) == 0:
            return added
        #print(f'saved = {self.saved}')
        #print(f'results = {results}')
        #assert 'displacements' not in results
        all_matched_results = self._get_matched_results(results)
        for result in all_matched_results:
            #self.log.warning(f'  add result = {result}')
            if result not in self.saved:
                self.saved.add(result)
                added.append(result)
        #print(f'saved = {self.saved}')
        #assert 'displacements' not in self.saved, self.saved
        return added

    def remove(self, results: str | list[str]) -> list[str]:
        """removes a list/str of results"""
        #self.log.warning(f'removing...results={results}')
        removed = []
        if len(results) == 0:
            return removed
        all_matched_results = self._get_matched_results(results)
        for result in all_matched_results:
            if result in self.saved:
                self.saved.remove(result)
                removed.append(result)
        #disable_set = set(results)
        #self.saved.difference(disable_set)
        return removed

    def _get_matched_results(self, results: str | list[str]) -> list[str]:
        """handles expansion of regexs"""
        if isinstance(results, str):
            results = [results]
        all_matched_results = []
        #assert 'displacements' not in all_matched_results
        #assert 'displacements' not in results
        for result in results:
            #print(result)
            if result in self.allowed:
                #self.log.warning(f'allowed - {result}')
                all_matched_results.append(result)
                continue

            resulti = _get_regex(result)
            #print('resulti =', resulti)
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
            #assert 'displacements' not in matched_results
        #assert 'displacements' not in all_matched_results
        return all_matched_results

    def _found_result(self, result: str) -> None:
        if result not in self.allowed:
            msg = f"result={result!r} is invalid; the name changed or it's a typo"
            raise RuntimeError(msg) # check line ~640 in op2_f06_common.py if this is a new result
        self.found.add(result)

    def update(self, results: list[str]) -> None:
        #assert 'displacements' not in results
        for result in results:
            self.add(result)

    #def add_found_result(self, result):
        #pass

    def __repr__(self) -> str:
        """defines the repr"""
        msg = 'ResultSet:\n'
        msg += ' results:\n'
        for result in sorted(self.allowed):
            if result in self.saved:
                msg += f'  {result}\n'
            else:
                msg += f'  {result} (disabled)\n'
        return msg

def _get_regex(result: str) -> str:
    r"""
    Tack on a word boundary if we have a * at the beginning of the regex

    Case Input           Output              Description
    ==== =============   =============       ============
    1    \b*strain       \b*strain           null because \b was found
    2    strain          \wstrain\w          add word boundaries (to only find strain exactly)
    3    *strain*        \w.*strain.*\w      any letter/number before/after the star
    4    *strain         \w.*strain          any letter/number before the star
    5    strain*         strain.*\w          any letter/number after the star

    Notes
    =====
    \w    letter/number [A-Za-z0-9_]
    \b    word boundary  -> trigger
    \.    period
    .     any character

    TODO: validate
    """
    #resulti = r'\w' + result if not result.startswith('*') else result  # old
    if '\b' in result:
        # Case 1
        # the user has gotten too fancy, so we'll let them do exactly what they want
        resulti = result
        #print(f'A: resulti = {resulti}')
        return resulti

    if '*' not in result:
        # Case 2 - add word boundaries to only find "result"
        resulti = fr'\w{result}\w'
        #print(f'B: resulti = {resulti!r}')
        return resulti

    # basically we replace * with .*
    # then we add a word boundary to either side
    wdot = r'\w.'  # \w or \.
    res_startswith_star = result.startswith('*')
    res_endswith_star = result.endswith('*')
    if res_startswith_star and res_endswith_star:
        #wdot = '\w'  # works
        resulti = f'{wdot}{result}'
        if result[-2] != '.':
            resulti = f'{resulti[:-1]}.*'
        #print(f'C: resulti = {resulti!r}')
        return resulti

    if res_startswith_star:
        resulti = f'{wdot}{result}'
        #print(f'D: resulti = {resulti!r}')
        return resulti

    # endswith
    resulti = result
    if result[-2] != '.':
        resulti = f'{result[:-1]}.*'
    #print(f'E: resulti = {resulti!r}')
    return resulti

def add_results_of_exact_type(all_results: list[str], result_word: str) -> list[str]:
    """
    Parameters
    ----------
    all_results : list[str]
        all the results to check against
    result_word : str
        'stress', 'strain', ...
        'stress' lets in stress.cquad4_stress, but not psd.cquad4_stress or grid_point_stress_discontinuities
        changed in v1.3.5; previously allowed things like displacement in (oops!) as well as the two mentioned cases

    Returns
    -------
    results_to_add : list[str]
        the results to add

    """
    results_to_add = []
    #_result_word = '_' + result_word
    result_word_dot = result_word + '.'
    for resulti in all_results:
        lower_resulti = resulti.lower()
        if result_word in lower_resulti and result_word_dot in lower_resulti:
            #print(result_word, resulti)
            # if _result_word not in lower_resulti:
            #     print('_word =', _result_word, lower_resulti)
            #flag = result_word in resulti.lower()
            #print(f'result_word={result_word!r} resulti={resulti}; flag={flag}')
            results_to_add.append(resulti)
    # stress_results = [result if 'stress' in result.lower() for result in all_results]
    return results_to_add
