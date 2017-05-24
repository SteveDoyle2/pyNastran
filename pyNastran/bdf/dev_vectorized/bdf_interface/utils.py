from six import string_types


def wipe_empty_fields(card):
    """
    Removes an trailing Nones from the card.
    Also converts empty strings to None.

    Parameters
    ----------
    card : list[str]
        the fields on the card as a list

    Returns
    -------
    short_card : List[str]
        the card with no trailing blank fields
    """
    cardB = []
    for field in card:
        if isinstance(field, string_types):
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