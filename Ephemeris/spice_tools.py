import spiceypy as spice
import numpy as np

def get_objects(filename, display = False):
    objects = spice.spkobj(filename)
    ids, names, tcs_sec, tcs_cal = [],[],[],[]
    n = 0

    if display:
        print('\nObjects in %s:' %filename)

    for o in objects:
        ids.append(o)
        tc_sec = spice.wnfetd(spice.spkcov(filename, ids[n]), n)
        tc_cal = [spice.timout(f, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB") for f in tc_sec]

        tcs_sec.append(tc_sec)
        tcs_cal.append(tc_cal)

        try:
            names.append(id2body(o))
        
        except:
            names.append('Unknown name')
        
        if display:
            print('id: %i\t\tname: %s\t\ttc: %s ---> %s' %(ids[-1], names[-1], tc_cal[0], tc_cal[1]))

    return ids, names, tcs_sec, tcs_cal

def id2body(id_):
    return spice.bodc2n(id_)

def tc2array(tcs, steps):
    arr = np.zeros((steps, 1))
    arr[:,0] = np.linspace(tcs[0], tcs[1], steps)
    return arr

def get_ephemeris_data(target, times, frame, observer):
    return np.array(spice.spkezr(target, times, frame, 'NONE', observer)[0])