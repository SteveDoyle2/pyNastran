"""
defines:
  - str_to_html(log_type, filename, lineno, msg)
"""
import datetime
import html

#message colors
DARK_ORANGE = '#EB9100'

COLORS = {
    'COMMAND' : 'blue',
    'ERROR' : 'Crimson',
    'DEBUG' : DARK_ORANGE,
    'WARNING' : 'purple',
    'INFO' : 'green',
}

def str_to_html(log_type, filename, lineno, msg):
    """
    Converts the message to html

    Parameters
    ----------
    color : str
        the HTML color
    log_type : str
        the message type
    filename : str
        the filename the message came from
    lineno : int
        the line number the message came from
    message : str
        the message

    Returns
    -------
    html_msg : str
        the HTML message

    """
    tim = datetime.datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')
    #print('log_type = %s' % log_type)
    #print('filename = %s' % filename)
    #print('lineno = %s' % lineno)
    #print('msg = %s' % msg)
    #assert isinstance(msg, str), msg
    msg = html.escape(msg)

    #message colors
    color = COLORS[log_type]

    if filename.endswith('.pyc'):
        filename = filename[:-1]
    html_msg = get_html_msg(color, tim, log_type, filename, lineno, msg)
    return html_msg

def get_html_msg(color, tim, log_type, filename, lineno, msg):
    """
    converts the message to html

    Parameters
    ----------
    color : str
        the HTML color
    time : str
        the time for the message
    log_type : str
        the message type
    filename : str
        the filename the message came from
    lineno : int
        the line number the message came from
    message : str
        the message

    Returns
    -------
    html_msg : str
        the HTML message
    """
    # log_type, filename, lineno, msg
    html_msg = r'<font color="%s"> %s %s : %s:%i</font> %s <br>' % (
        color, tim, log_type, filename, lineno, msg.replace('\n', '<br>'))
    return html_msg
