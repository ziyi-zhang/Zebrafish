%% fminsearch simplex
options = optimset('Display', 'iter');
[x, fval, exitflag, output] = fminsearch(f_realDot, [20, 20, 15, 4, 7], options);
imshowTiff(realDot); hold on; viscircles([x(2) x(1)], x(4));colormap pink;

%  Iteration   Func-count     min f(x)         Procedure
%      0            1      2.48175e+06         
%      1            6       2.4648e+06         initial simplex
%      2            7       2.4648e+06         reflect
%      3            9      1.66413e+06         expand
%      4           10      1.66413e+06         reflect
%      5           11      1.66413e+06         reflect
%      6           13           911584         expand
%      7           14           911584         reflect
%      8           16          40452.7         expand
%      9           17          40452.7         reflect
%     10           18          40452.7         reflect
%     11           19          40452.7         reflect
%     12           20          40452.7         reflect
%     13           22          40452.7         contract inside
%     14           24         -7147.56         reflect
%     15           25         -7147.56         reflect
%     16           27         -7147.56         contract outside
%     17           29         -10489.2         contract inside
%     18           31           -17391         contract inside
%     19           33         -30284.4         reflect
%     20           35         -34380.7         contract outside
%     21           37         -35299.5         reflect
%     22           39         -92476.3         reflect
%     23           40         -92476.3         reflect
%     24           42         -92476.3         contract inside
%     25           44         -92476.3         contract inside
%     26           46         -92476.3         contract inside
%     27           48          -113041         expand
%     28           50          -115447         reflect
%     29           52          -117597         reflect
%     30           53          -117597         reflect
%     31           54          -117597         reflect
%     32           56          -120824         reflect
%     33           58          -140449         reflect
%     34           60          -140449         contract inside
%     35           62          -140449         contract inside
%     36           63          -140449         reflect
%     37           65          -143907         contract inside
%     38           67          -147315         contract outside
%     39           69          -149174         contract inside
%     40           70          -149174         reflect
%     41           71          -149174         reflect
%     42           73          -150236         reflect
%     43           75          -151981         contract inside
%     44           77          -153227         contract inside
%     45           79          -153255         contract inside
%     46           80          -153255         reflect
%     47           82          -153255         contract inside
%     48           83          -153255         reflect
%     49           85          -153900         contract inside
%     50           87          -154445         contract inside
%     51           89          -154445         contract inside
%     52           91          -157453         reflect
%     53           93          -157453         contract inside
%     54           94          -157453         reflect
%     55           96          -157453         contract outside
%     56           97          -157453         reflect
%     57           98          -157453         reflect
%     58          100          -157453         contract inside
%     59          102          -157453         contract inside
%     60          104          -157453         contract inside
%     61          106          -157453         contract inside
%     62          108          -157453         contract inside
%     63          110          -157557         contract inside
%     64          111          -157557         reflect
%     65          113          -157904         reflect
%     66          115          -157925         reflect
%     67          117          -158016         reflect
%     68          119          -158016         contract inside
%     69          121          -158106         contract inside
%     70          123          -158150         reflect
%     71          125          -158150         contract inside
%     72          127          -158154         reflect
%     73          129          -158352         contract inside
%     74          130          -158352         reflect
%     75          131          -158352         reflect
%     76          133          -158573         expand
%     77          134          -158573         reflect
%     78          135          -158573         reflect
%     79          136          -158573         reflect
%     80          138          -158701         reflect
%     81          139          -158701         reflect
%     82          141          -158969         expand
%     83          142          -158969         reflect
%     84          143          -158969         reflect
%     85          144          -158969         reflect
%     86          146          -159469         expand
%     87          147          -159469         reflect
%     88          148          -159469         reflect
%     89          149          -159469         reflect
%     90          151          -159901         expand
%     91          153          -160137         expand
%     92          155          -160417         expand
%     93          156          -160417         reflect
%     94          158          -161653         expand
%     95          159          -161653         reflect
%     96          161          -163370         expand
%     97          162          -163370         reflect
%     98          164          -163370         contract inside
%     99          165          -163370         reflect
%    100          167          -163370         contract outside
%    101          168          -163370         reflect
%    102          169          -163370         reflect
%    103          170          -163370         reflect
%    104          172          -163370         contract inside
%    105          174          -163370         contract inside
%    106          175          -163370         reflect
%    107          176          -163370         reflect
%    108          178          -164145         expand
%    109          179          -164145         reflect
%    110          180          -164145         reflect
%    111          181          -164145         reflect
%    112          182          -164145         reflect
%    113          184          -164864         reflect
%    114          186          -164864         contract inside
%    115          187          -164864         reflect
%    116          188          -164864         reflect
%    117          190          -164864         contract inside
%    118          191          -164864         reflect
%    119          192          -164864         reflect
%    120          193          -164864         reflect
%    121          195          -165158         reflect
%    122          197          -165158         contract inside
%    123          198          -165158         reflect
%    124          200          -165158         contract inside
%    125          201          -165158         reflect
%    126          203          -165172         reflect
%    127          205          -165172         contract inside
%    128          207          -165520         expand
%    129          209          -165520         contract inside
%    130          211          -165520         contract inside
%    131          213          -165520         contract outside
%    132          214          -165520         reflect
%    133          216          -165520         contract inside
%    134          218          -165657         expand
%    135          219          -165657         reflect
%    136          221          -165659         reflect
%    137          223          -165919         reflect
%    138          225          -165919         contract inside
%    139          227          -166287         expand
%    140          229          -166287         contract inside
%    141          230          -166287         reflect
%    142          231          -166287         reflect
%    143          233          -166764         expand
%    144          235          -166764         contract inside
%    145          236          -166764         reflect
%    146          237          -166764         reflect
%    147          238          -166764         reflect
%    148          239          -166764         reflect
%    149          241          -167379         expand
%    150          242          -167379         reflect
%    151          243          -167379         reflect
%    152          244          -167379         reflect
%    153          246          -168110         expand
%    154          248          -168822         expand
%    155          249          -168822         reflect
%    156          250          -168822         reflect
%    157          251          -168822         reflect
%    158          252          -168822         reflect
%    159          253          -168822         reflect
%    160          255          -168822         contract inside
%    161          257          -169147         expand
%    162          259          -169266         reflect
%    163          261          -169266         contract inside
%    164          263          -169707         reflect
%    165          264          -169707         reflect
%    166          265          -169707         reflect
%    167          267          -169707         contract outside
%    168          269          -169707         contract outside
%    169          271          -169707         contract outside
%    170          273          -169707         contract inside
%    171          275          -169707         contract inside
%    172          276          -169707         reflect
%    173          278          -169707         contract inside
%    174          280          -169707         contract inside
%    175          282          -169731         contract inside
%    176          284          -169840         reflect
%    177          286          -169840         contract inside
%    178          287          -169840         reflect
%    179          288          -169840         reflect
%    180          289          -169840         reflect
%    181          291          -169840         contract inside
%    182          293          -169840         contract inside
%    183          295          -169840         contract inside
%    184          296          -169840         reflect
%    185          297          -169840         reflect
%    186          299          -169840         contract inside
%    187          300          -169840         reflect
%    188          301          -169840         reflect
%    189          303          -169840         contract inside
%    190          305          -169840         contract inside
%    191          306          -169840         reflect
%    192          308          -169840         contract inside
%    193          310          -169851         expand
%    194          311          -169851         reflect
%    195          313          -169851         contract inside
%    196          314          -169851         reflect
%    197          315          -169851         reflect
%    198          316          -169851         reflect
%    199          318          -169851         contract inside
%    200          320          -169851         contract inside
%    201          321          -169851         reflect
%    202          323          -169867         expand
%    203          324          -169867         reflect
%    204          325          -169867         reflect
%    205          326          -169867         reflect
%    206          328          -169874         reflect
%    207          330          -169893         expand
%    208          331          -169893         reflect
%    209          332          -169893         reflect
%    210          333          -169893         reflect
%    211          335          -169927         expand
%    212          336          -169927         reflect
%    213          338          -169963         expand
%    214          339          -169963         reflect
%    215          340          -169963         reflect
%    216          342          -170044         expand
%    217          343          -170044         reflect
%    218          345          -170045         reflect
%    219          347          -170061         reflect
%    220          349          -170061         contract inside
%    221          351          -170061         contract inside
%    222          353          -170061         contract inside
%    223          355          -170061         contract inside
%    224          357          -170067         reflect
%    225          359          -170082         reflect
%    226          361          -170082         contract inside
%    227          363          -170085         reflect
%    228          365          -170085         contract inside
%    229          367          -170085         contract inside
%    230          369          -170091         reflect
%    231          371          -170091         reflect
%    232          372          -170091         reflect
%    233          374          -170091         contract inside
%    234          376          -170091         contract inside
%    235          378          -170091         contract inside
%    236          380          -170094         reflect
%    237          382          -170109         expand
%    238          383          -170109         reflect
%    239          385          -170109         contract inside
%    240          387          -170109         contract inside
%    241          388          -170109         reflect
%    242          390          -170109         contract inside
%    243          392          -170109         contract inside
%    244          393          -170109         reflect
%    245          394          -170109         reflect
%    246          395          -170109         reflect
%    247          397          -170109         contract outside
%    248          398          -170109         reflect
%    249          400          -170109         contract inside
%    250          401          -170109         reflect
%    251          403          -170112         reflect
%    252          405          -170112         contract inside
%    253          407          -170112         contract inside
%    254          408          -170112         reflect
%    255          410          -170112         contract inside
%    256          412          -170112         contract inside
%    257          414          -170113         reflect
%    258          415          -170113         reflect
%    259          417          -170117         expand
%    260          419          -170117         contract inside
%    261          420          -170117         reflect
%    262          422          -170117         contract inside
%    263          424          -170117         contract outside
%    264          425          -170117         reflect
%    265          427          -170117         contract inside
%    266          429          -170117         contract inside
%    267          431          -170117         reflect
%    268          432          -170117         reflect
%    269          434          -170120         reflect
%    270          435          -170120         reflect
%    271          436          -170120         reflect
%    272          437          -170120         reflect
%    273          438          -170120         reflect
%    274          440          -170120         contract inside
%    275          441          -170120         reflect
%    276          443          -170120         contract inside
%    277          444          -170120         reflect
%    278          446          -170120         contract inside
%    279          448          -170120         contract inside
%    280          449          -170120         reflect
%    281          451          -170120         contract inside
%    282          453          -170120         contract inside
%    283          455          -170120         contract inside
%    284          457          -170120         contract inside
%    285          458          -170120         reflect
%    286          460          -170121         reflect
%    287          461          -170121         reflect
%    288          462          -170121         reflect
%    289          464          -170121         reflect
%    290          466          -170121         reflect
%    291          467          -170121         reflect
%    292          469          -170121         contract inside
%    293          471          -170121         contract inside
%    294          473          -170121         reflect
%    295          475          -170121         contract inside
%    296          477          -170121         contract inside
%    297          479          -170121         reflect
%    298          481          -170121         contract inside
%    299          482          -170121         reflect
%    300          483          -170121         reflect
%    301          485          -170121         contract inside
%    302          487          -170121         reflect
%    303          489          -170121         contract inside
%    304          491          -170121         reflect
%    305          493          -170122         reflect
%    306          495          -170122         contract outside
%    307          497          -170122         reflect
%    308          498          -170122         reflect
%    309          500          -170122         contract inside
%    310          502          -170122         contract inside
%    311          504          -170122         contract inside
%    312          506          -170122         contract outside
%    313          508          -170122         contract inside
%    314          510          -170122         reflect
%    315          512          -170122         expand
%    316          514          -170122         contract inside
%    317          516          -170122         contract inside
%    318          517          -170122         reflect
%    319          518          -170122         reflect
%    320          520          -170122         contract inside
%    321          521          -170122         reflect
%    322          523          -170122         contract inside
%    323          524          -170122         reflect
%    324          526          -170122         contract inside
%    325          528          -170122         contract inside
%    326          529          -170122         reflect
%    327          531          -170122         contract inside
%    328          532          -170122         reflect
%    329          534          -170122         contract inside
%    330          536          -170122         contract inside
%    331          538          -170122         contract inside
%    332          540          -170122         contract inside
%    333          541          -170122         reflect
%    334          543          -170122         reflect
%    335          545          -170122         contract inside
%    336          547          -170122         contract outside
%    337          549          -170122         contract inside
%    338          550          -170122         reflect
%    339          552          -170122         contract inside
%    340          554          -170122         contract inside
%    341          556          -170122         contract inside
%    342          557          -170122         reflect
%    343          559          -170122         contract outside
%    344          561          -170122         contract inside
%    345          563          -170122         contract inside
%    346          565          -170122         expand
%    347          566          -170122         reflect
%    348          568          -170122         contract inside
%    349          569          -170122         reflect
%    350          570          -170122         reflect
%    351          572          -170122         expand
%    352          573          -170122         reflect
%    353          575          -170122         contract inside
%    354          576          -170122         reflect
%    355          578          -170122         reflect
%    356          579          -170122         reflect
%    357          581          -170122         contract inside
%    358          583          -170122         contract inside
%    359          584          -170122         reflect
%    360          586          -170122         contract inside
%    361          587          -170122         reflect
%    362          589          -170122         reflect
%    363          591          -170122         contract inside
%    364          592          -170122         reflect
%    365          594          -170122         contract inside
%    366          596          -170122         contract outside
%    367          597          -170122         reflect
%    368          598          -170122         reflect
%    369          599          -170122         reflect
%    370          601          -170122         expand
%    371          603          -170122         contract inside
%    372          605          -170122         contract outside
%    373          607          -170122         contract inside
%    374          609          -170122         contract inside
%    375          611          -170122         expand
%    376          613          -170122         expand
%    377          615          -170122         contract inside
%    378          617          -170122         contract inside
%    379          618          -170122         reflect
%    380          620          -170122         reflect
%    381          622          -170122         expand
%    382          624          -170122         expand
%    383          625          -170122         reflect
%    384          626          -170122         reflect
%    385          627          -170122         reflect
%    386          629          -170122         expand
%    387          631          -170122         contract inside
%    388          632          -170122         reflect
%    389          633          -170122         reflect
%    390          635          -170122         expand
%    391          637          -170122         expand
%    392          639          -170122         expand
%    393          640          -170122         reflect
%    394          641          -170122         reflect
%    395          643          -170122         contract outside
%    396          645          -170122         contract inside
%    397          647          -170122         contract inside
%    398          649          -170122         contract inside
%    399          651          -170122         reflect
%    400          653          -170122         contract inside
%    401          655          -170122         contract inside
%    402          657          -170122         contract inside
%    403          659          -170122         contract inside
%    404          661          -170122         contract inside
%    405          662          -170122         reflect
%    406          663          -170122         reflect
%    407          665          -170122         reflect
%    408          666          -170122         reflect
%    409          667          -170122         reflect
%    410          669          -170122         contract inside
%    411          670          -170122         reflect
%    412          672          -170122         contract inside
%    413          674          -170122         contract inside
%    414          676          -170122         contract inside
%    415          678          -170122         expand
%    416          680          -170122         reflect
%    417          681          -170122         reflect
%    418          683          -170122         reflect
%    419          685          -170122         expand
%    420          687          -170122         reflect
%    421          689          -170122         contract inside
%    422          690          -170122         reflect
%    423          691          -170122         reflect
%    424          693          -170122         reflect
%    425          695          -170122         contract inside
%    426          697          -170122         expand
%    427          699          -170122         reflect
%    428          701          -170122         contract inside
%    429          702          -170122         reflect
%    430          704          -170122         reflect
%    431          706          -170122         reflect
%    432          707          -170122         reflect
%    433          709          -170122         contract outside
%    434          711          -170122         contract inside
%    435          713          -170122         contract inside
%    436          715          -170122         contract outside
%    437          717          -170122         contract inside
%    438          719          -170122         contract inside
%    439          720          -170122         reflect
%    440          722          -170122         contract inside
%    441          724          -170122         reflect
%    442          725          -170122         reflect
%    443          726          -170122         reflect
%    444          728          -170122         contract inside
%    445          730          -170122         reflect
%    446          732          -170122         contract outside
%    447          734          -170122         contract inside
%    448          735          -170122         reflect
%    449          737          -170122         reflect
%    450          739          -170122         contract inside
%    451          740          -170122         reflect
%    452          741          -170122         reflect
%    453          743          -170122         reflect
%    454          745          -170122         reflect
%    455          746          -170122         reflect
%    456          747          -170122         reflect
%    457          749          -170122         contract inside
%    458          751          -170122         contract outside
%    459          753          -170122         reflect
%    460          755          -170122         contract inside
%    461          757          -170122         reflect
%    462          759          -170122         contract inside
%    463          761          -170122         reflect
%    464          763          -170122         contract inside
%    465          765          -170122         contract inside
%    466          767          -170122         reflect
%    467          769          -170122         contract inside
%    468          771          -170122         contract inside
%    469          773          -170122         contract inside
 
% Optimization terminated:
%  the current x satisfies the termination criteria using OPTIONS.TolX of 1.000000e-04 
%  and F(X) satisfies the convergence criteria using OPTIONS.TolFun of 1.000000e-04 
% x =
% 
%   Columns 1 through 4
% 
%   16.315934428837963  16.442866125218949  17.573863631732436   4.234006052939736
% 
%   Column 5
% 
%    7.000000000596559

options = optimset('Display', 'iter');
[x, fval, exitflag, output] = fminsearch(f_testMat, [190, 310, 5, 18, 7], options);

%  Iteration   Func-count     min f(x)         Procedure
%      0            1           -74.25         
%      1            6         -161.143         initial simplex
%      2            8         -226.494         reflect
%      3           10         -268.851         expand
%      4           12         -308.087         reflect
%      5           13         -308.087         reflect
%      6           15         -308.087         contract outside
%      7           17         -308.087         contract inside
%      8           19         -308.087         contract inside
%      9           21         -313.224         contract inside
%     10           23         -323.769         reflect
%     11           24         -323.769         reflect
%     12           25         -323.769         reflect
%     13           27         -326.461         contract inside
%     14           29         -356.659         reflect
%     15           31         -356.659         contract inside
%     16           33         -356.659         contract inside
%     17           35         -356.659         contract outside
%     18           37         -356.659         contract inside
%     19           38         -356.659         reflect
%     20           39         -356.659         reflect
%     21           40         -356.659         reflect
%     22           42         -361.436         contract inside
%     23           43         -361.436         reflect
%     24           50         -361.436         shrink
%     25           52         -362.299         contract outside
%     26           53         -362.299         reflect
%     27           54         -362.299         reflect
%     28           56         -365.166         reflect
%     29           57         -365.166         reflect
%     30           59         -366.821         reflect
%     31           61         -366.821         contract inside
%     32           63         -366.821         contract inside
%     33           65         -366.821         contract inside
%     34           72         -366.821         shrink
%     35           74         -370.091         reflect
%     36           75         -370.091         reflect
%     37           77         -370.091         contract inside
%     38           79         -370.091         contract inside
%     39           81         -370.091         contract inside
%     40           82         -370.091         reflect
%     41           84         -370.091         contract inside
%     42           86         -370.091         contract inside
%     43           87         -370.091         reflect
%     44           89         -370.938         reflect
%     45           91         -370.938         contract inside
%     46           92         -370.938         reflect
%     47           94         -372.946         expand
%     48           95         -372.946         reflect
%     49           96         -372.946         reflect
%     50           98         -372.946         contract inside
%     51          100         -372.946         contract inside
%     52          102         -374.865         expand
%     53          104         -377.001         expand
%     54          106         -381.107         expand
%     55          107         -381.107         reflect
%     56          109         -385.176         expand
%     57          111         -387.544         expand
%     58          113          -389.48         expand
%     59          115         -390.777         reflect
%     60          117         -391.553         reflect
%     61          119         -391.553         contract inside
%     62          121         -391.553         contract inside
%     63          122         -391.553         reflect
%     64          124         -391.553         contract inside
%     65          126         -391.553         contract outside
%     66          128         -391.553         contract inside
%     67          129         -391.553         reflect
%     68          131         -391.897         reflect
%     69          132         -391.897         reflect
%     70          134         -393.336         reflect
%     71          135         -393.336         reflect
%     72          136         -393.336         reflect
%     73          138         -393.336         contract inside
%     74          139         -393.336         reflect
%     75          140         -393.336         reflect
%     76          141         -393.336         reflect
%     77          143         -393.544         expand
%     78          145         -393.544         contract outside
%     79          147         -393.544         contract outside
%     80          148         -393.544         reflect
%     81          150          -393.93         reflect
%     82          151          -393.93         reflect
%     83          153         -394.216         reflect
%     84          155         -394.269         contract inside
%     85          156         -394.269         reflect
%     86          157         -394.269         reflect
%     87          159         -394.337         expand
%     88          161         -394.614         reflect
%     89          162         -394.614         reflect
%     90          163         -394.614         reflect
%     91          164         -394.614         reflect
%     92          166         -394.769         reflect
%     93          168         -394.769         contract inside
%     94          170         -394.769         contract outside
%     95          172         -394.769         contract inside
%     96          173         -394.769         reflect
%     97          175         -394.804         reflect
%     98          177         -394.804         contract inside
%     99          179         -394.804         contract inside
%    100          180         -394.804         reflect
%    101          182         -394.843         expand
%    102          183         -394.843         reflect
%    103          185         -394.925         expand
%    104          186         -394.925         reflect
%    105          187         -394.925         reflect
%    106          189         -394.935         reflect
%    107          191         -394.993         reflect
%    108          193         -395.017         expand
%    109          194         -395.017         reflect
%    110          196         -395.017         contract inside
%    111          198         -395.088         reflect
%    112          199         -395.088         reflect
%    113          201         -395.088         contract inside
%    114          203         -395.088         contract outside
%    115          204         -395.088         reflect
%    116          206         -395.088         contract inside
%    117          208         -395.088         contract inside
%    118          210         -395.098         reflect
%    119          212         -395.111         reflect
%    120          213         -395.111         reflect
%    121          215         -395.165         reflect
%    122          217         -395.165         contract inside
%    123          219         -395.165         contract inside
%    124          221         -395.165         contract inside
%    125          222         -395.165         reflect
%    126          224         -395.165         contract inside
%    127          225         -395.165         reflect
%    128          227         -395.165         contract inside
%    129          229         -395.181         contract inside
%    130          230         -395.181         reflect
%    131          232         -395.181         contract outside
%    132          233         -395.181         reflect
%    133          234         -395.181         reflect
%    134          235         -395.181         reflect
%    135          237         -395.181         contract inside
%    136          238         -395.181         reflect
%    137          240         -395.185         reflect
%    138          242         -395.203         expand
%    139          244         -395.203         contract inside
%    140          246         -395.203         contract inside
%    141          248         -395.203         contract inside
%    142          250          -395.22         expand
%    143          251          -395.22         reflect
%    144          252          -395.22         reflect
%    145          254          -395.22         contract outside
%    146          256         -395.227         reflect
%    147          257         -395.227         reflect
%    148          259          -395.24         expand
%    149          260          -395.24         reflect
%    150          262         -395.256         expand
%    151          263         -395.256         reflect
%    152          265         -395.306         expand
%    153          266         -395.306         reflect
%    154          267         -395.306         reflect
%    155          269         -395.345         expand
%    156          271         -395.398         expand
%    157          273         -395.483         expand
%    158          274         -395.483         reflect
%    159          275         -395.483         reflect
%    160          277          -395.67         expand
%    161          278          -395.67         reflect
%    162          280         -395.819         expand
%    163          281         -395.819         reflect
%    164          282         -395.819         reflect
%    165          284         -396.312         expand
%    166          285         -396.312         reflect
%    167          286         -396.312         reflect
%    168          288         -396.842         expand
%    169          289         -396.842         reflect
%    170          290         -396.842         reflect
%    171          292         -396.963         expand
%    172          294         -397.918         expand
%    173          295         -397.918         reflect
%    174          297         -398.603         expand
%    175          298         -398.603         reflect
%    176          299         -398.603         reflect
%    177          301         -399.087         reflect
%    178          303         -399.187         reflect
%    179          305         -399.291         reflect
%    180          307         -399.291         contract inside
%    181          308         -399.291         reflect
%    182          310         -399.291         contract inside
%    183          312         -399.291         contract inside
%    184          314         -399.382         reflect
%    185          316         -399.382         contract inside
%    186          318         -399.382         contract outside
%    187          319         -399.382         reflect
%    188          321         -399.382         contract inside
%    189          323         -399.382         contract inside
%    190          325         -399.382         contract inside
%    191          327         -399.558         expand
%    192          329         -399.558         contract inside
%    193          331         -399.574         reflect
%    194          332         -399.574         reflect
%    195          334         -399.574         contract inside
%    196          335         -399.574         reflect
%    197          337         -399.574         contract inside
%    198          339         -399.574         contract inside
%    199          341         -399.574         contract inside
%    200          342         -399.574         reflect
%    201          344         -399.574         contract inside
%    202          346         -399.583         contract inside
%    203          348         -399.585         contract inside
%    204          350         -399.603         contract inside
%    205          352         -399.622         reflect
%    206          354         -399.622         contract inside
%    207          356         -399.622         contract inside
%    208          358         -399.622         contract inside
%    209          360         -399.622         contract inside
%    210          362         -399.622         contract inside
%    211          364         -399.622         contract outside
%    212          365         -399.622         reflect
%    213          367         -399.627         reflect
%    214          368         -399.627         reflect
%    215          370         -399.628         reflect
%    216          372         -399.628         contract inside
%    217          374         -399.637         reflect
%    218          376         -399.637         contract inside
%    219          378         -399.642         reflect
%    220          380         -399.642         contract inside
%    221          382         -399.642         contract inside
%    222          384         -399.642         contract inside
%    223          385         -399.642         reflect
%    224          386         -399.642         reflect
%    225          387         -399.642         reflect
%    226          389         -399.651         reflect
%    227          391         -399.651         contract inside
%    228          393         -399.651         contract inside
%    229          395         -399.651         contract inside
%    230          397         -399.651         contract inside
%    231          398         -399.651         reflect
%    232          399         -399.651         reflect
%    233          401         -399.651         contract inside
%    234          403         -399.651         contract outside
%    235          405         -399.651         contract inside
%    236          406         -399.651         reflect
%    237          408         -399.653         reflect
%    238          410         -399.658         expand
%    239          412         -399.658         contract inside
%    240          413         -399.658         reflect
%    241          415         -399.658         contract inside
%    242          417         -399.658         contract inside
%    243          418         -399.658         reflect
%    244          419         -399.658         reflect
%    245          421         -399.658         contract inside
%    246          423         -399.658         contract inside
%    247          425         -399.659         reflect
%    248          427         -399.659         contract inside
%    249          429         -399.659         contract inside
%    250          431         -399.659         contract inside
%    251          433         -399.659         reflect
%    252          434         -399.659         reflect
%    253          436         -399.659         contract inside
%    254          438         -399.659         contract inside
%    255          440          -399.66         reflect
%    256          442          -399.66         contract inside
%    257          444         -399.662         expand
%    258          446         -399.662         contract inside
%    259          447         -399.662         reflect
%    260          449         -399.662         contract inside
%    261          451         -399.662         contract inside
%    262          452         -399.662         reflect
%    263          454         -399.662         contract inside
%    264          455         -399.662         reflect
%    265          457         -399.662         contract inside
%    266          459         -399.662         contract inside
%    267          461         -399.662         contract inside
%    268          462         -399.662         reflect
%    269          464         -399.662         reflect
%    270          466         -399.662         contract inside
%    271          468         -399.662         contract inside
%    272          470         -399.662         reflect
%    273          471         -399.662         reflect
%    274          473         -399.662         contract inside
%    275          475         -399.662         contract inside
%    276          477         -399.662         contract inside
%    277          479         -399.662         contract inside
%    278          481         -399.662         contract inside
%    279          483         -399.662         reflect
%    280          484         -399.662         reflect
%    281          486         -399.662         contract inside
%    282          488         -399.662         contract outside
%    283          490         -399.662         contract inside
%    284          491         -399.662         reflect
%    285          492         -399.662         reflect
%    286          494         -399.662         contract inside
%    287          495         -399.662         reflect
%  
% Optimization terminated:
%  the current x satisfies the termination criteria using OPTIONS.TolX of 1.000000e-04 
%  and F(X) satisfies the convergence criteria using OPTIONS.TolFun of 1.000000e-04 
% 
% x
% 
% x =
% 
%    1.0e+02 *
% 
%   Columns 1 through 4
% 
%    2.001405668980957   3.000257495305228   0.047215957706743   0.198622768526376
% 
%   Column 5
% 
%    0.070000000985912

%% fminunc quasi-newton
options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
[x, fval, exitflag, output] = fminunc(f_testMat, [190, 310, 7, 17, 8], options);

%% particle swarm
options = optimoptions(@particleswarm, 'Display', 'iter');
[x, fval, exitflag, output] = particleswarm(f_realDot, 5, [1 1 1 2 1], [31 31 40 6 15], options);

%                                  Best            Mean     Stall
% Iteration     f-count            f(x)            f(x)    Iterations
%     0              50       2.299e+05             NaN        0
%     1             100       2.299e+05       2.147e+09        0
%     2             150       2.299e+05       6.791e+08        1
%     3             200       2.299e+05       1.264e+09        2
%     4             250       2.299e+05        9.55e+08        3
%     5             300       2.299e+05       9.767e+08        4
%     6             350       2.299e+05       8.059e+08        5
%     7             400       2.299e+05        8.27e+08        6
%     8             450       2.299e+05       2.271e+08        7
%     9             500      -2.284e+05       9.458e+05        0
%    10             550      -2.486e+05       1.799e+08        0
%    11             600      -2.486e+05       5.909e+07        1
%    12             650      -2.577e+05       1.201e+08        0
%    13             700       -3.27e+05       7.199e+05        0
%    14             750      -4.011e+05       6.341e+05        0
%    15             800      -4.048e+05        9.26e+05        0
%    16             850      -4.059e+05       8.856e+05        0
%    17             900      -4.059e+05       7.192e+05        1
%    18             950      -4.402e+05       5.073e+07        0
%    19            1000      -4.421e+05       1.139e+06        0
%    20            1050       -4.47e+05       5.576e+07        0
%    21            1100      -4.647e+05       6.707e+05        0
%    22            1150      -4.647e+05       5.723e+07        1
%    23            1200      -4.655e+05        1.08e+08        0
%    24            1250      -4.827e+05        8.97e+05        0
%    25            1300      -4.827e+05       5.979e+05        1
%    26            1350      -4.827e+05       6.557e+07        2
%    27            1400      -4.827e+05       7.597e+05        3
%    28            1450      -4.827e+05       2.083e+08        4
%    29            1500      -4.827e+05       1.201e+08        5
%    30            1550      -4.827e+05       7.001e+07        6
% 
%                                  Best            Mean     Stall
% Iteration     f-count            f(x)            f(x)    Iterations
%    31            1600      -4.827e+05       6.045e+05        7
%    32            1650      -4.827e+05       5.863e+05        8
%    33            1700      -4.839e+05       7.355e+05        0
%    34            1750      -4.846e+05       5.884e+05        0
%    35            1800      -4.852e+05       4.119e+05        0
%    36            1850      -4.875e+05       8.021e+05        0
%    37            1900      -4.908e+05       4.547e+05        0
%    38            1950      -4.975e+05       4.166e+05        0
%    39            2000      -4.975e+05       6.036e+07        1
%    40            2050      -4.975e+05       4.336e+05        2
%    41            2100      -4.975e+05       3.031e+05        0
%    42            2150      -4.983e+05         5.7e+05        0
%    43            2200      -4.984e+05       3.967e+05        0
%    44            2250      -4.987e+05       5.221e+05        0
%    45            2300      -4.999e+05       3.172e+05        0
%    46            2350      -4.999e+05        6.71e+05        1
%    47            2400      -4.999e+05       5.537e+05        2
%    48            2450      -4.999e+05       6.369e+07        3
%    49            2500      -5.005e+05       5.212e+05        0
%    50            2550      -5.005e+05       1.656e+08        1
%    51            2600      -5.005e+05       1.913e+05        2
%    52            2650      -5.005e+05        1.75e+05        3
%    53            2700      -5.005e+05       1.345e+08        4
%    54            2750      -5.005e+05       4.521e+05        5
%    55            2800      -5.005e+05       7.531e+05        6
%    56            2850      -5.005e+05       4.925e+05        7
%    57            2900      -5.007e+05       7.887e+05        0
%    58            2950      -5.007e+05       5.296e+07        1
%    59            3000      -5.008e+05       3.788e+05        0
%    60            3050      -5.008e+05       4.712e+05        0
% 
%                                  Best            Mean     Stall
% Iteration     f-count            f(x)            f(x)    Iterations
%    61            3100      -5.008e+05       7.565e+05        0
%    62            3150      -5.009e+05       6.697e+05        0
%    63            3200      -5.009e+05       5.934e+05        0
%    64            3250      -5.009e+05       9.313e+05        0
%    65            3300      -5.009e+05        7.99e+05        1
%    66            3350      -5.009e+05       4.112e+05        0
%    67            3400      -5.009e+05       5.431e+07        0
%    68            3450      -5.009e+05       4.915e+05        0
%    69            3500      -5.009e+05       3.998e+05        0
%    70            3550       -5.01e+05       6.724e+05        0
%    71            3600       -5.01e+05       3.542e+05        1
%    72            3650       -5.01e+05       5.417e+07        0
%    73            3700       -5.01e+05       5.706e+07        0
%    74            3750       -5.01e+05       5.068e+05        1
%    75            3800       -5.01e+05       3.194e+05        2
%    76            3850       -5.01e+05       5.839e+07        3
%    77            3900       -5.01e+05       1.135e+08        4
%    78            3950       -5.01e+05       6.551e+07        5
%    79            4000       -5.01e+05       4.117e+05        6
%    80            4050       -5.01e+05       6.004e+07        7
%    81            4100       -5.01e+05       5.002e+05        8
%    82            4150       -5.01e+05       4.167e+05        9
%    83            4200       -5.01e+05       3.658e+05        0
%    84            4250       -5.01e+05       5.475e+05        0
%    85            4300       -5.01e+05       3.923e+05        0
%    86            4350       -5.01e+05       7.145e+05        1
%    87            4400       -5.01e+05       6.901e+05        0
%    88            4450       -5.01e+05       6.263e+05        0
%    89            4500       -5.01e+05       4.252e+05        0
%    90            4550       -5.01e+05       5.256e+05        0
% 
%                                  Best            Mean     Stall
% Iteration     f-count            f(x)            f(x)    Iterations
%    91            4600       -5.01e+05       7.174e+05        0
%    92            4650       -5.01e+05        6.07e+05        0
%    93            4700       -5.01e+05       5.766e+05        0
%    94            4750       -5.01e+05       7.606e+05        1
%    95            4800       -5.01e+05       4.277e+05        2
%    96            4850       -5.01e+05       3.664e+05        3
%    97            4900       -5.01e+05       6.031e+07        0
%    98            4950       -5.01e+05       1.166e+08        0
%    99            5000       -5.01e+05       7.547e+05        1
%   100            5050       -5.01e+05       2.086e+05        2
%   101            5100       -5.01e+05       1.206e+05        3
%   102            5150       -5.01e+05       5.702e+07        4
%   103            5200       -5.01e+05       1.794e+08        5
%   104            5250       -5.01e+05       1.135e+08        6
%   105            5300       -5.01e+05       6.753e+07        7
%   106            5350       -5.01e+05       5.423e+07        0
%   107            5400       -5.01e+05       5.861e+07        0
%   108            5450       -5.01e+05       5.112e+05        0
%   109            5500       -5.01e+05       7.162e+05        0
%   110            5550       -5.01e+05       3.672e+05        0
%   111            5600       -5.01e+05       5.508e+05        0
%   112            5650       -5.01e+05       3.393e+05        0
% Optimization ended: relative change in the objective value 
% over the last OPTIONS.MaxStallIterations iterations is less than OPTIONS.FunctionTolerance.
% x
% 
% x =
% 
%   Columns 1 through 4
% 
%   16.482512180560107  16.515645057107992  18.374999650355743   3.815203604683403
% 
%   Column 5
% 
%    2.000000384008277
   
options = optimoptions(@particleswarm, 'Display', 'iter');
[x, fval, exitflag, output] = particleswarm(f_testMat, 5, [180 280 1 1 1], [220 320 40 40 40], options);

%                                  Best            Mean     Stall
% Iteration     f-count            f(x)            f(x)    Iterations
%     0              50             -85       9.019e+08        0
%     1             100             -85       9.607e+08        0
%     2             150             -85       4.668e+08        1
%     3             200            -133       1.342e+08        0
%     4             250            -167       1.342e+08        0
%     5             300            -167        8.59e+07        1
%     6             350            -242       1.315e+08        0
%     7             400            -242       4.295e+07        1
%     8             450            -270           93.59        0
%     9             500            -270       4.295e+07        1
%    10             550            -270       1.718e+08        2
%    11             600            -280           120.4        0
%    12             650            -280       8.948e+07        1
%    13             700            -280       4.295e+07        2
%    14             750            -287       4.474e+07        0
%    15             800            -287           149.9        1
%    16             850            -287       4.383e+07        2
%    17             900            -287           150.2        3
%    18             950            -287       4.474e+07        4
%    19            1000            -300           45.55        0
%    20            1050            -300           14.14        1
%    21            1100            -390           -80.6        0
%    22            1150            -390          -120.7        1
%    23            1200            -390          -201.2        2
%    24            1250            -438          -280.7        0
%    25            1300            -459          -341.7        0
%    26            1350            -479          -384.7        0
%    27            1400            -500          -413.6        0
%    28            1450            -503          -447.2        0
%    29            1500            -517            -468        0
%    30            1550            -520          -489.7        0
% 
%                                  Best            Mean     Stall
% Iteration     f-count            f(x)            f(x)    Iterations
%    31            1600            -520          -495.1        1
%    32            1650            -520          -489.1        2
%    33            1700            -520          -493.1        3
%    34            1750            -520          -494.1        4
%    35            1800            -520          -498.7        5
%    36            1850            -520          -493.6        6
%    37            1900            -520          -502.1        7
%    38            1950            -520          -506.6        8
%    39            2000            -520          -515.1        9
%    40            2050            -520          -515.1       10
%    41            2100            -520          -515.5       11
%    42            2150            -520          -518.6       12
%    43            2200            -520          -519.5       13
%    44            2250            -520          -518.9       14
%    45            2300            -520          -519.3       15
%    46            2350            -520          -518.8       16
%    47            2400            -520          -519.6       17
%    48            2450            -520          -519.7       18
%    49            2500            -520          -519.8       19
% Optimization ended: relative change in the objective value 
% over the last OPTIONS.MaxStallIterations iterations is less than OPTIONS.FunctionTolerance.
% x
% 
% x =
% 
%    1.0e+02 *
% 
%   Columns 1 through 4
% 
%    1.999349326574328   2.998607463343966   0.049762137282681   0.199603070415748
% 
%   Column 5
% 
%    0.101001958904583

%% Simulated Annealing
options = optimoptions(@simulannealbnd, 'Display', 'iter');
[x, fval, exitflag, output] = simulannealbnd(f_testMat, [190, 310, 5, 18, 7], [180 280 1 1 1], [220 320 40 40 40], options);

%                            Best        Current           Mean
% Iteration   f-count         f(x)         f(x)         temperature
%      0          1         -74.25         -74.25            100
%     10         11       -141.861       -141.861          56.88
%     20         21       -223.095       -223.095        34.0562
%     30         31       -223.095       -223.095        20.3907
%     40         41       -223.095       -219.437        12.2087
%     50         51       -223.095       -219.437        7.30977
%     60         61       -230.315       -230.315        4.37663
%     70         71       -347.052       -347.052        2.62045
%     80         81       -370.705       -370.705        1.56896
%     90         91       -435.895       -435.895       0.939395
%    100        101       -467.035       -467.035        0.56245
%    110        111       -475.349       -475.349        0.33676
%    120        121       -477.187       -477.187       0.201631
%    130        131       -480.741       -480.741       0.120724
%    140        141       -482.633       -482.633      0.0722817
%    150        151       -483.308       -483.308      0.0432777
%    160        161        -484.11        -484.11       0.025912
%    170        171       -485.359       -485.359      0.0155145
%    180        181       -485.815       -485.808     0.00928908
%    190        191       -486.353       -486.353     0.00556171
%    200        201       -486.687       -486.687        0.00333
%    210        211       -486.798       -486.798      0.0019938
%    220        221       -486.884       -486.884     0.00119376
%    230        231       -486.956       -486.956    0.000714748
%    240        241       -486.995       -486.995    0.000427946
%    250        251       -487.012       -487.012    0.000256227
%    260        261       -487.023       -487.023    0.000153413
%    270        271       -487.028       -487.028    9.18538e-05
%    280        281        -487.03        -487.03    5.49963e-05
% *  282        288       -487.031       -487.031        30.1632
%    290        296       -510.983       -510.983        20.0109
%    300        306       -510.983       -510.983        11.9813
% 
%                            Best        Current           Mean
% Iteration   f-count         f(x)         f(x)         temperature
%    310        316       -510.983       -510.983        7.17363
%    320        326       -510.983       -510.983        4.29512
%    330        336       -510.983       -510.983        2.57164
%    340        346       -510.983       -510.983        1.53974
%    350        356       -513.596       -513.596       0.921898
%    360        366       -519.739       -519.739       0.551975
%    370        376       -519.739       -519.739       0.330488
%    380        386       -519.985        -518.95       0.197875
%    390        396       -519.985       -519.945       0.118475
%    400        406           -520           -520      0.0709354
%    410        416           -520       -519.997      0.0424717
%    420        426           -520           -520      0.0254294
%    430        436           -520           -520      0.0152255
%    440        446           -520           -520     0.00911607
%    450        456           -520           -520     0.00545813
%    460        466           -520           -520     0.00326798
%    470        476           -520           -520     0.00195666
%    480        486           -520           -520     0.00117153
%    490        496           -520           -520    0.000701436
%    500        506           -520           -520    0.000419975
%    510        516           -520           -520    0.000251455
%    520        526           -520           -520    0.000150555
%    530        536           -520           -520     9.0143e-05
%    540        546           -520           -520    5.39719e-05
%    550        556           -520           -520     3.2315e-05
%    560        566           -520           -520    1.93482e-05
% *  561        572           -520           -520        44.3195