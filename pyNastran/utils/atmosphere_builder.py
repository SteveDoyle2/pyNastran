from typing import Any

import numpy as np
import matplotlib.pyplot as plt

from pyNastran.utils.atmosphere import (
    atm_equivalent_airspeed, atm_mach, atm_calibrated_airspeed,
    get_alt_for_eas_with_constant_mach,
    #get_alt_for_q_with_constant_mach,
    #get_alt_for_pressure,
    get_alt_for_mach_eas,
    get_mach_for_alt_eas,
    get_mach_for_alt_cas,
    #cas_to_mach,
    get_alt_for_cas_mach,
    #get_alt_for_mach_cas,
)


def build_envelope():
    alt_units = 'ft'
    eas_units = 'knots'
    cas_units = 'knots'

    xaxis = 'mach'
    #xaxis = 'cas'
    #xaxis = 'eas'
    #yaxis = 'alt'
    #yaxis = 'eas'
    yaxis = 'cas'
    data = [
        # const_type, (var_start, val_start), (var_end, val_end, name)),
        {
            'label': 'low_const_mach1',
            'tag': 'low_const_mach1',
            'const_type': 'eas',
            'start': {'alt': 0, 'mach': 0.2},
            'start_name': '0',

            'end':   {'alt': 40000,}, #[('alt', 40000), ('eas', '0')],
            'npoints': 50,
            'end_name': 'minQ_maxAlt',
        },
        {
            'label': 'const_alt2',
            'tag': 'const_alt2',
            'const_type': 'alt',
            'start': 'minQ_maxAlt',
            'end': {'alt': 40000, 'mach': 0.8},
            'npoints': 5,
            'end_name': 'maxQ_maxAlt',
        },
        #{
        #    'label': 'max_q1',
        #    'tag': 'max_q1',
        #    'const_type': 'eas',
        #    'start': {'alt': 0, 'eas': 225},
        #    'start_name': 'max_q0',
        #    'end': {'mach': 0.8,},
        #    'npoints': 10,
        #    'end_name': 'maxQ_maxMach',
        #},
        {
            'label': 'max_q2',
            'tag': 'max_q2',
            'const_type': 'eas',
            'start': {'alt': 0, 'eas': 250},
            'start_name': 'max_q0',
            #'end': {'alt': 40000,},
            'end': {'mach': 0.8,},
            'npoints': 10,
            'end_name': 'maxQ_minAlt',
        },
        {
            'label': 'min_alt',
            'tag': 'min_alt',
            'const_type': 'alt',
            'start': '0',
            'end': 'max_q0',
            'npoints': 10,
        },
        {
            'label': 'max_mach',
            'tag': 'max_mach',
            'const_type': 'mach',
            'start': 'maxQ_maxAlt',
            'end': 'maxQ_minAlt',
            'npoints': 10,
            'end_name': 'maxQ_maxAlt',
        },
        {
            'tag': 'max_cas_alt',
            'const_type': 'cas',
            'start': {'alt': 40000, 'mach': 0.8},
            'end': {'alt': 10000,},
            'npoints': 10,
        },
       #{
       #    'tag': 'max_cas_mach',
       #    'const_type': 'cas',
       #    'start': {'alt': 40000, 'mach': 0.8},
       #    'end': {'mach': 0.},
       #    'npoints': 10,
       #},
       #{
       #    'tag': 'cas_alt',
       #    'const_type': 'mach',
       #    'start': {'alt': 40000, 'cas': 200},
       #    'end': {'alt': 0.},
       #    'npoints': 10,
       #},
       #{
       #    'tag': 'cas_mach',
       #    'const_type': 'mach',
       #    'start': {'mach': 0.5, 'cas': 200},
       #    'end': {'alt': 0.},
       #    'npoints': 10,
       #},
    ]
    start = {}
    end = {}
    stored_data = {}
    mylines = []

    fig = plt.figure(1)
    ax = fig.gca()
    icolor = 0
    for idata, datai in enumerate(data):
        data_outi = _build_line(idata, datai,
                                stored_data,
                                alt_units, eas_units, cas_units)
        xdata = data_outi[xaxis]
        ydata = data_outi[yaxis]
        linestyle = datai.get('linestyle', '-')
        color = datai.get('color', f'C{idata}')
        label = datai.get('label', None)

        tag = datai['tag']
        alt = data_outi['alt']
        mach = data_outi['mach']
        eas = data_outi['eas']
        starti = {'alt': alt[0], 'mach': mach[0], 'eas': eas[0]}
        endi = {'alt': alt[-1], 'mach': mach[-1], 'eas': eas[-1]}

        assert tag not in start, tag
        start[tag] = starti
        end[tag] = endi
        ax.plot(xdata, ydata,
                linestyle=linestyle, label=label, color=color)
        mylines.append(data_outi)
        #print('-----------------')

    fend(ax, xaxis, yaxis)
    #-------------------------------------------------
    #print('---------------------------------------')
    scale = 1.1
    eas0 = start['low_const_mach1']['eas']
    data2 = [
        {
            'label': 'low_const_mach1-B',
            'tag': 'low_const_mach1-B',
            'const_type': 'eas',
            'start': {'alt': start['low_const_mach1']['alt'], 'eas': eas0/scale},
            #'start': {'mach': start['low_const_mach1']['mach'], 'eas': scale*eas0},
            'start_name': '0',

            'end':   {'mach': end['low_const_mach1']['mach'],},
            'npoints': 50,
            'end_name': 'minQ_maxAlt-B',
        },
        {
            'label': 'const_alt2-B',
            'tag': 'const_alt2-B',
            'const_type': 'alt',
            'start': 'minQ_maxAlt-B',
            'end': {'mach': 0.8*scale},
            'npoints': 5,
            'end_name': 'machMach_maxAlt-B',
        },
        {
            'label': 'max_q1-B',
            'tag': 'max_q1-B',
            'const_type': 'eas',
            'start': {'alt': 0, 'eas': 225*scale},
            'start_name': 'max_q0-B',
            'end': {'mach': 0.8*scale,},
            'npoints': 10,
            'end_name': 'maxQ_maxMach-B',
        },
    ]
    linestyle = '--'
    for idata, datai in enumerate(data2):
        data_outi = _build_line(idata, datai,
                                stored_data,
                                alt_units, eas_units, cas_units)
        xdata = data_outi[xaxis]
        ydata = data_outi[yaxis]
        #linestyle = datai.get('linestyle', '-')
        color = datai.get('color', f'C{idata}')
        label = datai.get('label', None)

        #print(color)
        alt = data_outi['alt']
        mach = data_outi['mach']
        eas = data_outi['eas']
        #print(f'**alt={alt}')
        #print(f'**eas={eas}')
        #starti = {'alt': alt[0], 'mach': mach[0], 'eas': eas[0]}
        #endi = {'alt': alt[-1], 'mach': mach[-1], 'eas': eas[-1]}
        #starts.append(starti)
        #ends.append(endi)
        ax.plot(xdata, ydata,
                linestyle=linestyle, label=label, color=color)
        mylines.append(data_outi)

    fend(ax, xaxis, yaxis)
    return mylines

def fend(ax, xaxis, yaxis):
    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)
    ax.grid(True)
    ax.legend()
    plt.show()


def _build_line(idata: int, data: dict[str, Any],
                stored_data,
                alt_units: str, eas_units: str, cas_units: str):
    #print(idata, data)
    tag = data['tag']
    const_type = data['const_type']
    start = data['start']
    end = data['end']

    npoints = data.get('npoints', 10)
    start_name = data.get('start_name', '')
    end_name = data.get('end_name', '')

    start_data = _point_to_data(tag, const_type, start, stored_data,
                                alt_units, eas_units, cas_units)
    #end_data = _point_to_data(const_type, end, stored_data)

    if not isinstance(start, str):
        assert len(start) == 2, (tag, start)

    if isinstance(end, str):
        end_data = _point_to_data(tag, const_type, end, stored_data,
                                  alt_units, eas_units, cas_units)
    else:
        assert len(end) >= 1, (tag, end)
        end_data = end

    #print(f'  start_data = {start_data}')
    if 'alt' in start_data and 'mach' in start_data:
        alt1, mach1 = start_data['alt'], start_data['mach']
    else:
        raise NotImplementedError((const_type, start_data))
    eas1 = atm_equivalent_airspeed(
        alt1, mach1, alt_units=alt_units, eas_units=eas_units)
    cas1 = atm_calibrated_airspeed(
        alt1, mach1, alt_units=alt_units, cas_units=cas_units)
    del start, end

    #print(f'  eas1 = {eas1:.0f} {eas_units}')
    if const_type == 'eas' and 'alt' in end_data:
        #print('  eas-alt', start_data, end_data)
        alt2 = end_data['alt']
        eas2 = eas1
        # get_alt_for_eas_with_constant_mach
        # get_alt_for_q_with_constant_mach
        # get_alt_for_pressure
        # get_alt_for_mach_eas(mach1, eas1, eas_units=eas_units)
        #mach2 = get_mach_for_alt_eas(
        #    alt2, eas2, alt_units=alt_units, eas_units=eas_units
        alt = np.linspace(alt1, alt2, num=npoints)
        mach = [get_mach_for_alt_eas(alti, eas1,
            alt_units=alt_units, eas_units=eas_units)
            for alti in alt]
    elif const_type == 'eas' and 'mach' in end_data:
        #print('  eas-mach', start_data, end_data)
        mach2 = end_data['mach']
        assert mach1 != mach2, (mach1, mach2)
        #print(f'eas-mach; mach2={mach2}') 
        eas2 = eas1
        mach = np.linspace(mach1, mach2, num=npoints)
        alt = [get_alt_for_eas_with_constant_mach(eas2, machi,
               alt_units=alt_units, velocity_units=eas_units)
               for machi in mach]
        #alt = np.array(alt).round(0)
        #print('alt =', alt)
        #print('mach =', mach)
        #print('eas =', eas)
    elif const_type == 'mach' and 'alt' in end_data:
        alt2 = end_data['alt']
        assert alt1 != alt2, (alt1, alt2)
        alt = np.linspace(alt1, alt2, num=npoints)
        mach = np.full(npoints, mach1, dtype='float64')

    elif const_type == 'alt' and 'mach' in end_data:
        mach2 = end_data['mach']
        mach = np.linspace(mach1, mach2, num=npoints)
        alt = np.full(npoints, alt1, dtype='float64')

    elif const_type == 'cas' and 'alt' in end_data:
        alt2 = end_data['alt']
        cas2 = cas1
        alt = np.linspace(alt1, alt2, num=npoints)
        mach = [get_alt_for_cas_mach(cas2, alti, alt_units=alt_units,
                                     cas_units=cas_units)
               for alti in alt]
    elif const_type == 'cas' and 'mach' in end_data:
        mach2 = end_data['mach']
        mach = np.linspace(mach1, mach2, num=npoints)
        cas2 = cas1
        alt = [get_alt_for_cas_mach(
                  cas2, machi, alt_units=alt_units,
                  cas_units=cas_units) for machi in mach]
    #elif const_type == 'mach' and 'alt' in end_data:
    else:
        raise NotImplementedError((const_type, end_data))

    eas = [atm_equivalent_airspeed(alti, machi, alt_units=alt_units, eas_units=eas_units)
           for alti, machi in zip(alt, mach)]


    cas = [atm_calibrated_airspeed(alti, machi, alt_units=alt_units, cas_units=cas_units)
           for alti, machi in zip(alt, mach)]
    
    alt = np.array(alt) #.round(0)
    mach = np.array(mach) #.round(3)
    eas = np.array(eas) #.round(0)
    cas = np.array(cas) #.round(0)
    #print('alt =', alt)
    #print('mach =', mach)
    #print('eas =', eas)

    alt2 = alt[-1]
    mach2 = mach[-1]
    eas2 = eas[-1]
    cas2 = cas[-1]

    if 0:  # pragma: no cover
        alt1 = float(round(alt1, 0))
        mach1 = float(round(mach1, 3))
        eas1 = float(round(eas1, 0))
        cas1 = float(round(cas1, 0))
        
        alt2 = float(round(alt2, 0))
        mach2 = float(round(mach2, 3))
        eas2 = float(round(eas2, 0))
        cas2 = float(round(cas2, 0))

    start_data = {
        'alt': alt1, 'mach': mach1,
        'eas': eas1, 'cas': cas1, }
    if start_name:
        #print(f'  saving start_name {start_name!r}; {start_data}')
        stored_data[start_name] = start_data

    end_data = {
        'alt': alt2, 'mach': mach2,
        'eas': eas2, 'cas': cas2,}
    if end_name:
        #print(f'  saving end_name {end_name!r}; {end_data}')
        stored_data[end_name] = end_data

    #const_type = data['const_type']
    data_out = {
        'alt': alt,
        'mach': mach,
        'eas': eas,
        'cas': cas,
    }
    return data_out

def _alt_mach(alt: float, mach: float, **kwargs):
    return alt, mach

def _alt_eas(alt: float, eas: float, **kwargs):
    alt_units = kwargs['alt_units']
    eas_units = kwargs['eas_units']
    mach = get_mach_for_alt_eas(
        alt, eas, alt_units=alt_units, eas_units=eas_units)
    return alt, mach

def _alt_cas(alt: float, cas: float, **kwargs):
    alt_units = kwargs['alt_units']
    cas_units = kwargs['cas_units']
    mach = get_mach_for_alt_cas(
        alt, cas, alt_units=alt_units, cas_units=cas_units)
    return alt, mach

def _eas_mach(eas: float, mach: float, **kwargs):
    alt_units = kwargs['alt_units']
    eas_units = kwargs['eas_units']
    alt = get_alt_for_eas_with_constant_mach(
        eas, mach, alt_units=alt_units, velocity_units=eas_units)
    return alt, mach

def _cas_mach(cas: float, mach: float, **kwargs):
    alt_units = kwargs['alt_units']
    cas_units = kwargs['cas_units']
    alt = get_alt_for_cas_mach(
        cas, mach, alt_units=alt_units, cas_units=cas_units)
    return alt, mach

def _cas_eas(cas: float, mach: float, **kwargs):
    cas_units = kwargs['cas_units']
    eas_units = kwargs['eas_units']
    #get_alt_mach_for_cas_eas(cas, eas, cas_units=cas_units, eas_units=eas_units)
    raise NotImplementedError('cas_eas')
    return alt, mach

def _point_to_data(tag: str,
                   const_type: str,
                   point_data,
                   stored_data: dict[str, tuple[float, float]],
                   alt_units: str, eas_units: str,
                   cas_units: str) -> tuple[float, float]:
    """
    maps 2 arbitrary quantities (alt, mach, eas, cas) to
    (alt, mach)
    """
    if isinstance(point_data, str):
        outi = stored_data[point_data]
        assert isinstance(outi, dict), outi
        assert len(outi) >= 2, (tag, outi)
        return outi

    assert isinstance(point_data, dict), point_data
    assert len(point_data) == 2, point_data
    
    arg_type = []
    value = []
    for arg_typei, valuei in sorted(point_data.items()):
        arg_type.append(arg_typei)
        value.append(valuei)
    arg_type1, arg_type2 = arg_type[:2]
    val1, val2 = value[:2]
    assert arg_type1 != arg_type2

    # alt, mach -> (eas, cas)
    func_map = {
        ('alt', 'mach'): _alt_mach,
        ('alt', 'eas'): _alt_eas,
        ('alt', 'cas'): _alt_cas,

        ('cas', 'eas'): _cas_eas,
        ('eas', 'mach'): _eas_mach,
        ('cas', 'mach'): _cas_mach,
        #(const_type, 'alt', 'mach'): passer,
        #{'eas': atm_equivalent_airspeed, atm_calibrated_airspeed)},
    }
    func = func_map[(arg_type1, arg_type2)]
    alt, mach = func(
        val1, val2,
        alt_units=alt_units,
        eas_units=eas_units,
        cas_units=cas_units)

    eas = atm_equivalent_airspeed(
        alt, mach, alt_units=alt_units, eas_units=eas_units)
    cas = atm_calibrated_airspeed(
        alt, mach, alt_units=alt_units, cas_units=cas_units)

    #if 0:  # pragma: no cover
    #    alt = round(alt, 0)
    #    mach = round(mach, 3)
    #    eas = round(eas, 1)
    #    cas =round(cas, 1)

    out = {
        'alt': float(alt),
        'mach': float(mach),
        'eas': float(eas),
        'cas': float(cas),
    }
    return out # alt, mach
    


if __name__ == '__main__':
    build_envelope()
