BDF Headers
===========

Version
-------
You've probably seen the following in your BDF:
```
$ pyNastran: version: msc
```

That's a shorthand for:
```
model = BDF(mode='msc')
```

What's neat the flag in the deck takes priority, so the user can set their version without modifying the code and it will work properly.

Allowable versions include:

 - msc (default)
 - nx
 - optistruct
 - mystran
 - zona

Punch
-----

Like version:
```
$ pyNastran: punch=True
```
overwrites the flag specified in the BDF() tag (or more likely adds it because it was left empty).

Encoding
---------

Encoding is an flag that's useful for decoding difficult characters.  latin1, cp1252, and utf8 are very common encodings.
```
$ pyNastran: encoding=cp1252
```


Skip Cards
----------
skip_cards is incredibly powerful when you have duplicate ids or just want to limit cards in the GUI.
```
$ pyNastran: skip_cards = CBUSH, PBUSH
```
This will delete the elements in the group from the file.  This is a general method, so:
```
$ pyNastran: skip elements=1:10 15 16
$ pyNastran: skip properties=1:10 15 16
$ pyNastran: skip materials=1:10 15 16
...
```
work too.

Skip
----
Skip is more intended for use in the GUI, but is also nice as a way to simplify the model.
```
$ pyNastran: skip elements=1:10 15 16
```
This will delete the elements in the group from the file.  This is a general method, so:
```
$ pyNastran: skip elements=1:10 15 16
$ pyNastran: skip properties=1:10 15 16
$ pyNastran: skip materials=1:10 15 16
...
```
work too.
