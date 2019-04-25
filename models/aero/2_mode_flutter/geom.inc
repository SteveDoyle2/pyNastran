PARAM    POST    -2
param,postext,yes
PARAM   PRTMAXIM YES
EIGRL    1                       20      0
$ Direct Text Input for Bulk Data
$ Elements and Element Properties for region : plate
PSHELL   1       1      .01      1               1
$ Pset: "plate" will be imported as: "pshell.1"
CQUAD4   1       1       1       2       31      30
CQUAD4   2       1       2       3       32      31
CQUAD4   3       1       3       4       33      32
CQUAD4   4       1       4       5       34      33
CQUAD4   5       1       5       6       35      34
CQUAD4   6       1       6       7       36      35
CQUAD4   7       1       7       8       37      36
CQUAD4   8       1       8       9       38      37
CQUAD4   9       1       9       10      39      38
CQUAD4   10      1       10      11      40      39
CQUAD4   11      1       11      12      41      40
CQUAD4   12      1       12      13      42      41
CQUAD4   13      1       13      14      43      42
CQUAD4   14      1       14      15      44      43
CQUAD4   15      1       15      16      45      44
CQUAD4   16      1       16      17      46      45
CQUAD4   17      1       17      18      47      46
CQUAD4   18      1       18      19      48      47
CQUAD4   19      1       19      20      49      48
CQUAD4   20      1       20      21      50      49
CQUAD4   21      1       21      22      51      50
CQUAD4   22      1       22      23      52      51
CQUAD4   23      1       23      24      53      52
CQUAD4   24      1       24      25      54      53
CQUAD4   25      1       25      26      55      54
CQUAD4   26      1       26      27      56      55
CQUAD4   27      1       27      28      57      56
CQUAD4   28      1       28      29      58      57
CQUAD4   29      1       30      31      60      59
CQUAD4   30      1       31      32      61      60
CQUAD4   31      1       32      33      62      61
CQUAD4   32      1       33      34      63      62
CQUAD4   33      1       34      35      64      63
CQUAD4   34      1       35      36      65      64
CQUAD4   35      1       36      37      66      65
CQUAD4   36      1       37      38      67      66
CQUAD4   37      1       38      39      68      67
CQUAD4   38      1       39      40      69      68
CQUAD4   39      1       40      41      70      69
CQUAD4   40      1       41      42      71      70
CQUAD4   41      1       42      43      72      71
CQUAD4   42      1       43      44      73      72
CQUAD4   43      1       44      45      74      73
CQUAD4   44      1       45      46      75      74
CQUAD4   45      1       46      47      76      75
CQUAD4   46      1       47      48      77      76
CQUAD4   47      1       48      49      78      77
CQUAD4   48      1       49      50      79      78
CQUAD4   49      1       50      51      80      79
CQUAD4   50      1       51      52      81      80
CQUAD4   51      1       52      53      82      81
CQUAD4   52      1       53      54      83      82
CQUAD4   53      1       54      55      84      83
CQUAD4   54      1       55      56      85      84
CQUAD4   55      1       56      57      86      85
CQUAD4   56      1       57      58      87      86
CQUAD4   57      1       59      60      89      88
CQUAD4   58      1       60      61      90      89
CQUAD4   59      1       61      62      91      90
CQUAD4   60      1       62      63      92      91
CQUAD4   61      1       63      64      93      92
CQUAD4   62      1       64      65      94      93
CQUAD4   63      1       65      66      95      94
CQUAD4   64      1       66      67      96      95
CQUAD4   65      1       67      68      97      96
CQUAD4   66      1       68      69      98      97
CQUAD4   67      1       69      70      99      98
CQUAD4   68      1       70      71      100     99
CQUAD4   69      1       71      72      101     100
CQUAD4   70      1       72      73      102     101
CQUAD4   71      1       73      74      103     102
CQUAD4   72      1       74      75      104     103
CQUAD4   73      1       75      76      105     104
CQUAD4   74      1       76      77      106     105
CQUAD4   75      1       77      78      107     106
CQUAD4   76      1       78      79      108     107
CQUAD4   77      1       79      80      109     108
CQUAD4   78      1       80      81      110     109
CQUAD4   79      1       81      82      111     110
CQUAD4   80      1       82      83      112     111
CQUAD4   81      1       83      84      113     112
CQUAD4   82      1       84      85      114     113
CQUAD4   83      1       85      86      115     114
CQUAD4   84      1       86      87      116     115
$ Referenced Material Records
$ Material Record : Aluminum
$ Description of Material : Date: 29-Nov-17           Time: 14:51:22
$ E = 6.89E10 Pa
$ G = ???
$ nu = 0.3
$ rho = 2800 kg/m^3?
MAT1     1      6.89+10         .3      2800.
$ Nodes of the Entire Model
GRID     1               0.      0.      0.
GRID     2               0.     .357     0.
GRID     3               0.     .714     0.
GRID     4               0.     1.071    0.
GRID     5               0.     1.429    0.
GRID     6               0.     1.786    0.
GRID     7               0.     2.143    0.
GRID     8               0.     2.5      0.
GRID     9               0.     2.857    0.
GRID     10              0.     3.214    0.
GRID     11              0.     3.571    0.
GRID     12              0.     3.929    0.
GRID     13              0.     4.286    0.
GRID     14              0.     4.643    0.
GRID     15              0.     5.       0.
GRID     16              0.     5.357    0.
GRID     17              0.     5.714    0.
GRID     18              0.     6.071    0.
GRID     19              0.     6.429    0.
GRID     20              0.     6.786    0.
GRID     21              0.     7.143    0.
GRID     22              0.     7.5      0.
GRID     23              0.     7.857    0.
GRID     24              0.     8.214    0.
GRID     25              0.     8.571    0.
GRID     26              0.     8.929    0.
GRID     27              0.     9.286    0.
GRID     28              0.     9.643    0.
GRID     29              0.     10.      0.
GRID     30             .333     0.      0.
GRID     31             .333    .357     0.
GRID     32             .333    .714     0.
GRID     33             .333    1.071    0.
GRID     34             .333    1.429    0.
GRID     35             .333    1.786    0.
GRID     36             .333    2.143    0.
GRID     37             .333    2.5      0.
GRID     38             .333    2.857    0.
GRID     39             .333    3.214    0.
GRID     40             .333    3.571    0.
GRID     41             .333    3.929    0.
GRID     42             .333    4.286    0.
GRID     43             .333    4.643    0.
GRID     44             .333    5.       0.
GRID     45             .333    5.357    0.
GRID     46             .333    5.714    0.
GRID     47             .333    6.071    0.
GRID     48             .333    6.429    0.
GRID     49             .333    6.786    0.
GRID     50             .333    7.143    0.
GRID     51             .333    7.5      0.
GRID     52             .333    7.857    0.
GRID     53             .333    8.214    0.
GRID     54             .333    8.571    0.
GRID     55             .333    8.929    0.
GRID     56             .333    9.286    0.
GRID     57             .333    9.643    0.
GRID     58             .333    10.      0.
GRID     59             .667     0.      0.
GRID     60             .667    .357     0.
GRID     61             .667    .714     0.
GRID     62             .667    1.071    0.
GRID     63             .667    1.429    0.
GRID     64             .667    1.786    0.
GRID     65             .667    2.143    0.
GRID     66             .667    2.5      0.
GRID     67             .667    2.857    0.
GRID     68             .667    3.214    0.
GRID     69             .667    3.571    0.
GRID     70             .667    3.929    0.
GRID     71             .667    4.286    0.
GRID     72             .667    4.643    0.
GRID     73             .667    5.       0.
GRID     74             .667    5.357    0.
GRID     75             .667    5.714    0.
GRID     76             .667    6.071    0.
GRID     77             .667    6.429    0.
GRID     78             .667    6.786    0.
GRID     79             .667    7.143    0.
GRID     80             .667    7.5      0.
GRID     81             .667    7.857    0.
GRID     82             .667    8.214    0.
GRID     83             .667    8.571    0.
GRID     84             .667    8.929    0.
GRID     85             .667    9.286    0.
GRID     86             .667    9.643    0.
GRID     87             .667    10.      0.
GRID     88             1.       0.      0.
GRID     89             1.      .357     0.
GRID     90             1.      .714     0.
GRID     91             1.      1.071    0.
GRID     92             1.      1.429    0.
GRID     93             1.      1.786    0.
GRID     94             1.      2.143    0.
GRID     95             1.      2.5      0.
GRID     96             1.      2.857    0.
GRID     97             1.      3.214    0.
GRID     98             1.      3.571    0.
GRID     99             1.      3.929    0.
GRID     100            1.      4.286    0.
GRID     101            1.      4.643    0.
GRID     102            1.      5.       0.
GRID     103            1.      5.357    0.
GRID     104            1.      5.714    0.
GRID     105            1.      6.071    0.
GRID     106            1.      6.429    0.
GRID     107            1.      6.786    0.
GRID     108            1.      7.143    0.
GRID     109            1.      7.5      0.
GRID     110            1.      7.857    0.
GRID     111            1.      8.214    0.
GRID     112            1.      8.571    0.
GRID     113            1.      8.929    0.
GRID     114            1.      9.286    0.
GRID     115            1.      9.643    0.
GRID     116            1.      10.      0.
$ Loads for Load Case : Default
$ Referenced Coordinate Frames

