import sys
import platform

#if platform.system=='Windows':


def get_graphic_card_properties():
    import dbus
    bus = dbus.SystemBus()
    hal_manager_object = bus.get_object(
        'org.freedesktop.Hal', '/org/freedesktop/Hal/Manager')
    prop = 'pci.device_class'
    for device in hal_manager_object.get_dbus_method('GetAllDevices', 'org.freedesktop.Hal.Manager')():
        dev = bus.get_object('org.freedesktop.Hal', device)
        interface = dbus.Interface(
            dev, dbus_interface='org.freedesktop.Hal.Device')
        if interface.PropertyExists(prop):
            if interface.GetProperty(prop) == 3:
                # we return the properties of the first device in the list
                # with a pci.device_class == 3 (should check if several such devs...
                return interface.GetAllProperties()

if 0:
    dic = get_graphic_card_properties()
    for key, value in dic.iteritems():
        print("%s : %s" % (key, value))


def makeLog():
    msg = ''
    msg += "sys.version           = %s\n" % (sys.version)
    msg += "sys.version_info      = %s\n" % (str(sys.version_info))
    msg += "machine               = %s\n" % (platform.machine())
    msg += "platform              = %s\n" % (platform.platform())
    msg += "processor             = %s\n" % (platform.processor())
    msg += "architecure           = %s\n" % (str(platform.architecture()))
    #msg += "os          = %s\n" %(platform.os())
    msg += "python_branch         = %s\n" % (platform.python_branch())
    msg += "python_revision       = %s\n" % (platform.python_revision())
    msg += "win32_ver             = %s\n" % (str(platform.win32_ver()))
    msg += "version               = %s\n" % (platform.version())
    msg += "uname                 = %s\n" % (str(platform.uname()))
    msg += "system                = %s\n" % (platform.system())
    msg += "python_build          = %s\n" % (str(platform.python_build()))
    msg += "python_compiler       = %s\n" % (platform.python_compiler())
    msg += "python_implementation = %s\n" % (platform.python_implementation())
    msg += "system                = %s\n" % (platform.system())
    #msg += "system_alias          = %s\n" %(platform.system_alias())
    msg += "mac_ver               = %s\n" % (str(platform.mac_ver()))
    msg += "linux_distribution    = %s\n" % (
        str(platform.linux_distribution()))
    msg += "libc_ver              = %s\n" % (str(platform.libc_ver()))
    print msg
    f = open('pyNastran.log', 'w')
    f.write(msg)
    f.close()

if __name__=='__main__':
    makeLog()
