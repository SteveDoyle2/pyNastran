from typing import TextIO
import numpy as np
from cpylog import SimpleLogger
from pyNastran.utils import object_attributes, object_methods


class CampbellData:
    def __init__(self, solution: int,
                 cddata_list: list[dict[int, np.ndarray]]):
        #print(f'solution = {solution:d}')
        self.solution = solution
        self.cddata_list: list[dict[int, np.ndarray]] = cddata_list

    def __eq__(self, obj) -> bool:
        return True

    def is_sort1(self) -> bool:
        return True

    def is_sort2(self) -> bool:
        return False

    def is_real(self) -> bool:
        return False

    def is_complex(self) -> bool:
        return True

    def export_to_hdf5(self, result_group, log: SimpleLogger):
        pass

    def write_f06(self, f06_file: TextIO):
        msg = ''
        for isub, resii in enumerate(self.cddata_list):
            for ii, resiii in resii.items():
                msg += f'{isub},{ii}: {resiii.tolist()}\n'
        f06_file.write(msg)
        # print(msg)

    def get_stats(self, short: bool=False) -> list[str]:
        msg = [
            '  type=%s solution=%s\n' % (self.__class__.__name__, self.solution),
            '  cddata_list\n'
        ]
        #self.plot()
        return msg

    def object_attributes(self, mode: str='public', keys_to_skip=None,
                          filter_properties: bool=False):
        """
        List the names of attributes of a class as strings. Returns public
        attributes as default.

        Parameters
        ----------
        obj : instance
            the object for checking
        mode : str
            defines what kind of attributes will be listed
            * 'public' - names that do not begin with underscore
            * 'private' - names that begin with single underscore
            * 'both' - private and public
            * 'all' - all attributes that are defined for the object
        keys_to_skip : list[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        attribute_names : list[str]
            sorted list of the names of attributes of a given type or None
            if the mode is wrong
        """
        return object_attributes(self, mode=mode, keys_to_skip=keys_to_skip,
                                 filter_properties=filter_properties)

    def object_methods(self, mode: str='public', keys_to_skip=None):
        """
        List the names of methods of a class as strings. Returns public methods
        as default.

        Parameters
        ----------
        obj : instance
            the object for checking
        mode : str
            defines what kind of methods will be listed
            * "public" - names that do not begin with underscore
            * "private" - names that begin with single underscore
            * "both" - private and public
            * "all" - all methods that are defined for the object
        keys_to_skip : list[str]; default=None -> []
            names to not consider to avoid deprecation warnings

        Returns
        -------
        method : list[str]
            sorted list of the names of methods of a given type
            or None if the mode is wrong
        """
        return object_methods(self, mode=mode, keys_to_skip=keys_to_skip)

    def plot(self, ifig: int=1, show: bool=True) -> None:
        import matplotlib.pyplot as plt
        for data in self.cddata_list:
            # 1: 'RPM',
            # 2: 'eigenfreq',
            # 3: 'Lehr',
            # 4: 'real_eig',
            # 5: 'imag_eig',
            # 6: 'whirl_dir',
            # 7: 'converted_freq',
            # 8: 'whirl_code',
            rpm = data['RPM']
            eigenfreq = data['eigenfreq']
            Lehr = data['Lehr']
            # ['RPM', 'eigenfreq', 'Lehr',
            #  'real_eig', 'whirl_dir', 'converted_freq',
            #  'whirl_code']
            converted_freq = data['converted_freq']
            plt.figure(ifig)
            plt.plot(rpm, eigenfreq)  # RPM vs. eigenfreq
            # plt.grid(True)

            plt.figure(ifig+1)
            plt.plot(rpm, Lehr)  # RPM vs. Lehr
            plt.grid(True)

            real_eig = data['real_eig']
            plt.figure(ifig+2)
            plt.plot(rpm, real_eig)  # RPM vs. real_eig
            plt.grid(True)

            # imag_eig = data['imag_eig']
            # plt.figure(ifig+3)
            # plt.plot(rpm, imag_eig)  # RPM vs. imag_eig
            # plt.grid(True)

            plt.figure(ifig+4)
            plt.plot(rpm, converted_freq)  # RPM vs. converted_freq
            plt.grid(True)
        if show:  # pragma: no cover
            plt.show()
