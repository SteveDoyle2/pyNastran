from pyNastran.utils import is_string


def wipe_empty_fields(card):
    """
    Removes an trailing Nones from the card.
    Also converts empty strings to None.

    :param card:        the fields on the card as a list
    :returns shortCard: the card with no trailing blank fields
    """
    cardB = []
    for field in card:
        if is_string(field):
            field = field.strip()
            if field == '':
                field = None
        cardB.append(field)

    i = 0
    iMax = 0
    while i < len(card):
        if cardB[i] is not None:
            iMax = i
        i += 1
    return cardB[:iMax + 1]