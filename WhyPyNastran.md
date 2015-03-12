# Why I develop the Software #
I (Steve Doyle) am an aerospace engineer for a small aerospace company ([M4 Engineering](http://www.m4-engineering.com/)) based out of Long Beach, CA.  I do a lot of aerodynamic analysis, but also work with Nastran every now and then.  We often do multi-disciplinary analysis and many of the tools we use on a daily basis are developed in-house.

The problem with some of these internal tools is it's difficult to improve a tool with so much legacy.  It's not practical to regularly rewrite your widely used libraries, so the software becomes hard to maintain.  Additionally, a fully capable Nastran reader requires way more funding than any company should have to invest.  I eventually realized the the only way to actually get a fully capable tool, I'd have to use an open-source tool.

Since the ones out there (e.g. [pynastran](http://www.koders.com/info.aspx?c=ProjectInfo&pid=9MZG4YNSDBPLZY5Q11KE8N11SF&s=cdef%3aparser) and [FeResPost](http://www.ferespost.be/)) were not what I was looking for (despite doing some very impressive things with their code), I decided to open-source part of my master's work (which is why there's also a Cart3D viewer).  I would have stopped there, but I realized in order for the project to grow, I needed to get some users and eventually developers.  I added a few features and when I realized I enjoyed doing it, I spent a lot of time making it easy to add new features.

At this point (at least in the trunk version) there are ~210 cards and many of the OP2 tables are supported.  I also was able to get Al Danial to write an ASCII/binary OP4 reader, who is also more knowledgeable about the inner workings of Nastran than I am.  The software isn't perfect, but despite the excessive amount of Nastran tests there are (e.g. the tpl directory), it's actually hard to find errors.  Additionally, there are a few people who make tickets for the harder to catch bugs.

The one catch about having a capable Nastran reader is the API will change slightly over time.  I try to avoid it, but it will happen and if something needs to be scrapped because it's buggy, it will be (e.g. the old F06 Reader).

So in summary:

# Why use pyNastran #
  * it's object oriented (this makes using it much easier)
  * the goal is to read everything accurately, which requires a lot of extra work up front, but it drastically reduces errors
  * it's written in [Python](http://www.python.org/) and works in Python 2 and Python 3, so no expensive licenses are required (e.g. Matlab)
  * it's free
  * by being open-source it can be improved by people who a minor update to fit their specific use case
  * you can debug you own code and fix it.  That said, if you want to release it as part of your own software, you need to release the fix (it's part of the license)
  * you can make a [ticket](http://code.google.com/p/pynastran/issues/list) and it will get fixed (unlike the tickets I make :)
  * it's really easy to add additional features (once you know your way around)
  * it hooks in well with other programs (e.g. the stress/deflection margin calculator)
  * it will work in much larger programs (e.g. [NASA's OpenMDAO](http://www.openmdao.org)) and even supports the OpenMDAO [optimization format](http://code.google.com/p/pynastran/wiki/Optimization) that's used in the their NastranWrapper
  * it supports large/small/CSV/tab formatted cards
  * it cards write out in small field format with minimal loss of precision


So ask questions, post in the forum, make tickets.  It's not going to do what you want it to unless you ask.  If the wiki isn't clear (or right), just let me know.

# Alternatives #
  * [IMAT](http://www.ata-e.com/downloads/imat/) is a very capable MATLAB toolbox developed by [ATA-E Engineering](http://www.ata-e.com/).
  * [NastranWrapper](https://github.com/OpenMDAO-Plugins/nastranwrapper) developed by NASA as part of [OpenMDAO](http://www.openmdao.org).  It's open-source, and meant to be used as part of OpenMDAO, which is also open-source.  My understanding is that it's meant as a quick find-replace tool to be be used as part of optimization, whereas pyNastran's goal is to read everything.  It's also very well done.
  * [FeResPost](http://www.ferespost.be/) a GPL Ruby-based extension (Windows/Linux) or a COM component or .NET assembly (Windows).  I think this works with VBA.  It's also in French.  I haven't used it, but it actually has real documentation and at a glance seems pretty capable.
  * [PyNASTRAN.py](http://translate.google.com/translate?hl=en&sl=de&u=http://mbi-wiki.uni-wuppertal.de/2010/07/&prev=/search%3Fq%3DPYNASTRAN%26start%3D40%26hl%3Den%26client%3Dfirefox-beta%26hs%3DOh8%26sa%3DN%26rls%3Dorg.mozilla:en-US:official%26biw%3D1079%26bih%3D588%26prmd%3Dimvns&sa=X&ei=abjqT6Vt8YDYBePu8bQB&ved=0CFQQ7gEwATgo) is the original pyNastran!  I accidentally stole the name :)  It's got a basic reader and viewer.  I was never able to get it to work.

If you have any others, let me know and it will be added.