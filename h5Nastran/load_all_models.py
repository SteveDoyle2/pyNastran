import os
import os.path

from collections import defaultdict

from h5Nastran import H5Nastran, pyNastranReadBdfError, pyNastranWriteBdfError

# this file verifies that h5Nastran can read all the bdf files included with pyNastran
# any bdf cards that are not supported are printed out


# path to pyNastran models, you will need to point to the correct path here
models_path = r'P:\redmond\pyNastran\models'

bdfs = []

for dirpath, dirnames, filenames in os.walk(models_path):
    for filename in [f for f in filenames if f.endswith('.bdf') or f.endswith('.dat')]:
        bdfs.append(os.path.join(dirpath, filename))

unsupported_cards = defaultdict(int)

for bdf in bdfs:
    db = H5Nastran('dummy.h5', 'w', in_memory=True)
    try:
        db.load_bdf(bdf)
    except (pyNastranReadBdfError, pyNastranWriteBdfError):
        db.close()
        continue

    for card_id in db._unsupported_bdf_cards:
        unsupported_cards[card_id] += 1

    db.close()

unsupported_cards = [(card_id, unsupported_cards[card_id]) for card_id in unsupported_cards.keys()]
unsupported_cards = sorted(unsupported_cards, key=lambda x: x[1])

for card in unsupported_cards:
    print(card)

print(len(unsupported_cards))



