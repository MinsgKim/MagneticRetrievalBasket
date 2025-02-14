#ifndef RTW_HEADER_MagneticBasket_Simscape_Optimizer_h_
#define RTW_HEADER_MagneticBasket_Simscape_Optimizer_h_
#ifndef MagneticBasket_Simscape_Optimizer_COMMON_INCLUDES_
#define MagneticBasket_Simscape_Optimizer_COMMON_INCLUDES_
#include <stdlib.h>
#include "sl_AsyncioQueue/AsyncioQueueCAPI.h"
#include "rtwtypes.h"
#include "sigstream_rtw.h"
#include "simtarget/slSimTgtSigstreamRTW.h"
#include "simtarget/slSimTgtSlioCoreRTW.h"
#include "simtarget/slSimTgtSlioClientsRTW.h"
#include "simtarget/slSimTgtSlioSdiRTW.h"
#include "simstruc.h"
#include "fixedpoint.h"
#include "raccel.h"
#include "slsv_diagnostic_codegen_c_api.h"
#include "rt_logging_simtarget.h"
#include "dt_info.h"
#include "ext_work.h"
#include "nesl_rtw.h"
#include "MagneticBasket_Simscape_Optimizer_dda62cd9_1_gateway.h"
#endif
#include "MagneticBasket_Simscape_Optimizer_types.h"
#include <stddef.h>
#include "rtw_modelmap_simtarget.h"
#include "rt_defines.h"
#include <string.h>
#include "rtGetInf.h"
#include "rt_nonfinite.h"
#define MODEL_NAME MagneticBasket_Simscape_Optimizer
#define NSAMPLE_TIMES (2) 
#define NINPUTS (0)       
#define NOUTPUTS (0)     
#define NBLOCKIO (333) 
#define NUM_ZC_EVENTS (0) 
#ifndef NCSTATES
#define NCSTATES (18)   
#elif NCSTATES != 18
#error Invalid specification of NCSTATES defined in compiler command
#endif
#ifndef rtmGetDataMapInfo
#define rtmGetDataMapInfo(rtm) (*rt_dataMapInfoPtr)
#endif
#ifndef rtmSetDataMapInfo
#define rtmSetDataMapInfo(rtm, val) (rt_dataMapInfoPtr = &val)
#endif
#ifndef IN_RACCEL_MAIN
#endif
typedef struct { real_T aay54cmgia [ 18 ] ; real_T nfyrk0gene [ 154 ] ;
real_T gdleddk2bu [ 7 ] ; real_T o3cnonqf2y [ 7 ] ; real_T bjyrxtyjst [ 3 ] ;
real_T ecnsyysan1 [ 3 ] ; real_T d11srgkxv1 [ 3 ] ; real_T oywwaak55z ;
real_T n3uyejqyiv [ 3 ] ; real_T dg2tz20qh3 [ 9 ] ; real_T jd0epngxde ;
real_T i1wwep0quu ; real_T jytdis5y4s [ 4 ] ; real_T dvlkqx1wrk [ 3 ] ;
real_T axlew0gacz [ 3 ] ; real_T ffk0uso1jv [ 9 ] ; real_T j2caeelflf [ 9 ] ;
real_T gzdxak0rgb [ 9 ] ; real_T nhpvoyt4rs [ 3 ] ; real_T mnl5wynqkw ;
real_T n5zmerun4k [ 3 ] ; real_T i210rcby1e [ 6 ] ; real_T njlotg3eyk [ 3 ] ;
real_T ciemy5f0t3 [ 3 ] ; real_T i3fsmpf0l4 [ 3 ] ; real_T p1t00nt2c4 ;
real_T lkk1p3hkce ; real_T j0nzmiwajs ; real_T p4b33uxmxl [ 4 ] ; real_T
kerf1yzlz0 [ 3 ] ; real_T dxk5dr2guu [ 3 ] ; real_T gfwkvsxkxh [ 3 ] ; real_T
jmgra3bzio ; real_T n43m55vtpu ; real_T dm0nsvocoi [ 3 ] ; real_T pzoqxgrplf
[ 3 ] ; real_T pi3p3kites [ 3 ] ; real_T j4tq1cqvap [ 3 ] ; real_T d4sl1kvckp
; real_T k0e5tb5t4r [ 3 ] ; real_T mrse3gkxz4 [ 3 ] ; real_T hjkncigkim [ 3 ]
; real_T bv20qr1sq1 [ 3 ] ; real_T jknqpcyts5 ; real_T fo44k43jik [ 3 ] ;
real_T leoubbnbua [ 9 ] ; real_T eb5u3a2x55 ; real_T eded4se0e2 ; real_T
epqt4glrc2 [ 4 ] ; real_T huyyr23tvs [ 3 ] ; real_T a5iwna5ofw [ 3 ] ; real_T
n5qkyzyc0p [ 9 ] ; real_T cviiadazje [ 9 ] ; real_T hrdkh50xl2 [ 9 ] ; real_T
p4bi154apu [ 3 ] ; real_T ollqza00rh ; real_T kypr3ufomq [ 3 ] ; real_T
ljklcjruel [ 6 ] ; real_T giyavk0qxh [ 3 ] ; real_T iwp0nxhxfg [ 3 ] ; real_T
cwwkrda1qf [ 3 ] ; real_T ojzjawyp2p ; real_T j22ln0sz4g ; real_T ltesgv5ex5
; real_T e4k3bfcn1z [ 4 ] ; real_T fi20a5zfv5 [ 3 ] ; real_T oexwmjcbqo [ 3 ]
; real_T ljdiomnr0e [ 3 ] ; real_T nj4zxwlnaz ; real_T b4eypdlrbg ; real_T
mswvoafvyk [ 3 ] ; real_T ga3a4xdbbc [ 3 ] ; real_T eecvbm2tyj [ 3 ] ; real_T
bkxjvl2wwl ; real_T jcdegmljub [ 3 ] ; real_T dakqbtmvna [ 3 ] ; real_T
j54ujtxdmj [ 3 ] ; real_T ol0mge5uht [ 3 ] ; real_T fusdzrgc4y [ 3 ] ; real_T
n2cw123fsj ; real_T gtc4x0vmrj [ 3 ] ; real_T igdgu1b2hk [ 9 ] ; real_T
lntlwwhs4g ; real_T lvldddhyjb ; real_T j1izfy2apg [ 4 ] ; real_T f2e0m2twoo
[ 3 ] ; real_T pctih4a2gv [ 3 ] ; real_T d15mjte2mm [ 9 ] ; real_T o0sl3javac
[ 9 ] ; real_T lk2g50sct1 [ 9 ] ; real_T kucavif4r4 [ 3 ] ; real_T fd02dku0jj
; real_T j3daqzj5si [ 3 ] ; real_T e30x0ooogi [ 6 ] ; real_T dgv155wtql [ 3 ]
; real_T bcueeadxdi [ 3 ] ; real_T nalsmrc4um [ 3 ] ; real_T gmj04jp0xd ;
real_T ij2rjc05nz ; real_T pyxwiqcuay [ 4 ] ; real_T ndcqoyp0vx [ 3 ] ;
real_T hv30h5tsxp ; real_T iqlb144ibi ; real_T pl32irduyc [ 3 ] ; real_T
lgifpnih0e [ 3 ] ; real_T l0eiwonrmi [ 3 ] ; real_T hrtkouuc3i [ 3 ] ; real_T
gmo4itlisv ; real_T nexxpkogi1 [ 3 ] ; real_T idwieepnv4 ; real_T h2vbgw30tf
[ 3 ] ; real_T dndiwk2b1s [ 3 ] ; real_T gvavv3jbs0 [ 3 ] ; real_T dw12uizmd2
[ 3 ] ; real_T bdswxd3b20 [ 3 ] ; real_T fh04srz5g0 ; real_T mkzi43umox [ 3 ]
; real_T cpzvyphms1 [ 9 ] ; real_T oc1jrzdnhe ; real_T hlnjhscgik ; real_T
fkzq14ohjb [ 4 ] ; real_T cafqksoz00 [ 3 ] ; real_T c402o0md5n [ 3 ] ; real_T
hwppgfjqwy [ 9 ] ; real_T i1lnvldtbv [ 9 ] ; real_T ahw1dydhgo [ 9 ] ; real_T
ort4y5zlaf [ 3 ] ; real_T f1arg5nwa3 ; real_T eckq3h2tzb [ 3 ] ; real_T
bsaqv0u04c [ 6 ] ; real_T alzuotumkg [ 3 ] ; real_T ley50ghvfj [ 3 ] ; real_T
b0dwpzqyrm [ 3 ] ; real_T fbdf2tllgu ; real_T my44vuapi0 ; real_T k3faf5atq0
[ 4 ] ; real_T cg3eiq4i4s [ 3 ] ; real_T encqrxuu3y ; real_T bdozk1fpny ;
real_T nbrrv2d3zg [ 3 ] ; real_T cuxk4j3qys [ 3 ] ; real_T eegu2ofo43 [ 3 ] ;
real_T hacw1qj3aj [ 3 ] ; real_T pxg5b2slaw ; real_T h4h1ivsmcr [ 3 ] ;
real_T jqzltcnumt ; real_T etawd0whn4 [ 3 ] ; real_T cbmxpbd2iu [ 3 ] ;
real_T l1wqaovs4e [ 3 ] ; real_T f34tk41pk1 [ 3 ] ; real_T kup2makpml [ 3 ] ;
real_T d4it13dar5 ; real_T dl4tl2eodk [ 3 ] ; real_T obk5xjqd1g [ 9 ] ;
real_T p00fanezjt ; real_T b4b10et3s3 ; real_T flzta24bjr [ 4 ] ; real_T
dv2svmenmy [ 3 ] ; real_T odlamrju54 [ 3 ] ; real_T p3ru4cddgi [ 9 ] ; real_T
bqjv3qidt0 [ 9 ] ; real_T lkoo4fzueb [ 9 ] ; real_T gxha2ko2vz [ 3 ] ; real_T
lffwcdbemm ; real_T chig3jpwbn [ 3 ] ; real_T fcsl1vvaio [ 6 ] ; real_T
pu15hoqemb [ 3 ] ; real_T emdpdwkzg3 [ 3 ] ; real_T lxihdocxga [ 3 ] ; real_T
nfak5a1vn4 ; real_T ntt11gc034 ; real_T hvd2sf2hag [ 4 ] ; real_T bpewjnssvt
[ 3 ] ; real_T j1ydxe2usi ; real_T ceb1rzrjkf ; real_T lyhkcvfx5d [ 3 ] ;
real_T enoycngmy3 [ 3 ] ; real_T mw4hhpohry [ 3 ] ; real_T bhsuwz2d2l [ 3 ] ;
real_T ezvzwoavxd ; real_T msqcnwfh2s [ 3 ] ; real_T mp054zeo2w ; real_T
pfosbwya4i [ 3 ] ; real_T l01qgkxm4q [ 3 ] ; real_T babxk0cwue [ 3 ] ; real_T
fzbgmnan4p [ 3 ] ; real_T oflmpiaez1 [ 3 ] ; real_T efmfmwthqv ; real_T
ivaqfjshfj [ 3 ] ; real_T fzqka3pwrj [ 9 ] ; real_T egvr4ozujh ; real_T
cjrpkdpszw ; real_T logpc2hca5 [ 4 ] ; real_T ksrswfmi3x [ 3 ] ; real_T
d2uj0bf0jq [ 3 ] ; real_T joimmmlrrj [ 9 ] ; real_T aaw4gtn4yc [ 9 ] ; real_T
jnz14mmcr3 [ 9 ] ; real_T le0lvwzojd [ 3 ] ; real_T k3vxq1m344 ; real_T
e3hmodlujb [ 3 ] ; real_T b3yujzbge2 [ 6 ] ; real_T gf5n1xgexk [ 3 ] ; real_T
e3buvmif1c [ 3 ] ; real_T j4fldanieq [ 3 ] ; real_T cdgxtvlobx ; real_T
lbmtuwi01h ; real_T ieyrepmb05 ; real_T lyjnnzstpa [ 4 ] ; real_T favu2vq25u
[ 3 ] ; real_T h1zdr1aani [ 3 ] ; real_T m4bi5zptkh [ 3 ] ; real_T cab5adcq1r
; real_T avcgwoxwhd ; real_T f4nnctaymu [ 3 ] ; real_T mcasp4wuog [ 3 ] ;
real_T mrdivqgmvm [ 3 ] ; real_T g3vpjafupq [ 3 ] ; real_T btylsma43j ;
real_T my424itvki [ 3 ] ; real_T fmqgelgm0i [ 3 ] ; real_T gmivsxqevf [ 3 ] ;
real_T cxwqk2amda [ 3 ] ; real_T jgkbpstejl ; real_T evftq3sk0q [ 3 ] ;
real_T hnzvyzmoy4 [ 9 ] ; real_T gtu1mufolx ; real_T ifl2ongx1d ; real_T
a3qlbfbxxo [ 4 ] ; real_T glnk1qsnt1 [ 3 ] ; real_T li0pbutf3d [ 3 ] ; real_T
nzexd0tfyn [ 9 ] ; real_T lrkvtt341w [ 9 ] ; real_T kptzapguuy [ 9 ] ; real_T
nwco0rdsmq [ 3 ] ; real_T j2h4dmgm1o ; real_T kdzk5onxbh [ 3 ] ; real_T
oiy31a51av [ 6 ] ; real_T athb42mz2d [ 3 ] ; real_T exqp2uohe5 [ 3 ] ; real_T
o2znllyvoq [ 3 ] ; real_T bbtvepo4p4 ; real_T kq5inrtfnq ; real_T c00wpmk13d
; real_T huacygyj0o [ 4 ] ; real_T gg2oxfw35t [ 3 ] ; real_T lperucgciu [ 3 ]
; real_T jqxy4sgz1r [ 3 ] ; real_T aizglk3yjb ; real_T avjauevoe4 ; real_T
azxdea4e00 [ 3 ] ; real_T afxh4mt3gu [ 3 ] ; real_T gdhtjuaiqp [ 3 ] ; real_T
e0gxuicxzo ; real_T kbih5351xy [ 3 ] ; real_T j2tnneegmr [ 3 ] ; real_T
c513o5pg5r [ 4 ] ; real_T pfmi5qvzhf [ 4 ] ; real_T gi35bo1nw4 [ 4 ] ; real_T
iii5xyl5wn [ 4 ] ; real_T mjqwg1vruy [ 4 ] ; real_T eo0yhyzcnn [ 4 ] ; real_T
c5dvbg4axm [ 4 ] ; real_T bxbu5xamaf [ 4 ] ; real_T foqtyz1zpt [ 4 ] ; real_T
gy3u0julyl [ 4 ] ; real_T hytxiurw0s [ 4 ] ; real_T nubcraskrw [ 4 ] ; real_T
lp0cxari03 [ 4 ] ; real_T phx5lsta5w [ 4 ] ; real_T hoap0zr3eu [ 4 ] ; real_T
enp5f3s002 [ 4 ] ; real_T nsl31ekm0f [ 4 ] ; real_T hrqvwki2jl [ 4 ] ; real_T
dypcby02yn [ 4 ] ; real_T a100glwqkz [ 4 ] ; real_T c3twv1nzdj [ 4 ] ; real_T
iurh4h1cla [ 4 ] ; real_T har0a2wrmi [ 4 ] ; real_T ib2ieuf2an [ 4 ] ; real_T
ojpovrodvr [ 4 ] ; real_T fcgmpow5uz [ 4 ] ; real_T jylntbxfpj [ 4 ] ; real_T
g1yo2p3pxt [ 4 ] ; real_T az4eo3zh5o [ 4 ] ; real_T ch1tbnhyon [ 4 ] ; real_T
eoseoj2f3b [ 4 ] ; real_T aigpe1li2p [ 4 ] ; real_T pbn2nzk2k3 [ 4 ] ; real_T
hffhwq1imf [ 4 ] ; real_T gbxuu4wtte [ 4 ] ; real_T k5nlme3kns [ 4 ] ; real_T
j52cfj2d1w [ 4 ] ; real_T ftvz2ekipg [ 4 ] ; real_T awgskjjmsb [ 4 ] ; real_T
fl4rvtkarn [ 4 ] ; real_T jp22vgjwwm [ 4 ] ; real_T pzpjv5s5n2 [ 4 ] ; real_T
nsmqj1yv3h [ 3 ] ; real_T mdrxplc5qx [ 3 ] ; real_T idcqc0uibn [ 3 ] ; real_T
lkwtwvkutl [ 3 ] ; real_T dunjuwexfg [ 3 ] ; real_T oimhqzeh2t [ 3 ] ; real_T
bcudxny0c0 [ 3 ] ; real_T ibott2adsg [ 3 ] ; real_T ovchpiv2ws [ 3 ] ; real_T
eejq0hbv30 [ 3 ] ; real_T joumqp0zpt [ 3 ] ; real_T l02fn5ndrj [ 3 ] ; real_T
emllridil0 [ 3 ] ; real_T d1zqfcz4zo [ 3 ] ; real_T bjhxaihh05 [ 3 ] ; real_T
jltrdwyzmw [ 3 ] ; real_T fc15cwanga [ 3 ] ; real_T lebemiabnp [ 3 ] ; real_T
p5ihfgrspp [ 3 ] ; real_T itx5qkpva4 [ 3 ] ; real_T e4tk4xaa2e [ 3 ] ; real_T
mdxnlbsfvm [ 3 ] ; real_T fvtfdqea2r [ 3 ] ; real_T inhs5zftab [ 3 ] ; real_T
hrudfjbz3i [ 3 ] ; real_T cschdiniut [ 3 ] ; real_T o4wni3515w [ 3 ] ; real_T
p5hpxiz5u3 [ 3 ] ; real_T h2ica3jpzg [ 3 ] ; real_T jqlckttaam [ 3 ] ; real_T
npob1uuqas [ 3 ] ; real_T itxt4ytwah [ 3 ] ; real_T nl3crkyc5i [ 3 ] ; real_T
kdrism1ogn [ 3 ] ; real_T nuca2ykly0 [ 3 ] ; } B ; typedef struct { real_T
ik1trwr3gq [ 2 ] ; real_T m0rr22lxax [ 2 ] ; real_T jhm0cqabg3 [ 2 ] ; real_T
oq35dv23vw [ 2 ] ; real_T lwhmmfv2ps [ 2 ] ; real_T iyivz4czyo [ 2 ] ; real_T
ls5rwnwukt [ 2 ] ; real_T mebb30zlfb [ 2 ] ; real_T kmluc3bnxm [ 2 ] ; real_T
ardpmlcavq [ 2 ] ; real_T o2te45t5iz [ 2 ] ; real_T ceunniigow [ 2 ] ; real_T
by4chuqxsd [ 2 ] ; real_T pjygrikwyo [ 2 ] ; real_T dbasyf4pfi [ 2 ] ; real_T
nj2tur0kwz [ 2 ] ; real_T f1cbvtmv22 [ 2 ] ; real_T mh3gjthtbl [ 2 ] ; real_T
jwgslfbtic [ 2 ] ; real_T ezstggpf2j [ 2 ] ; real_T pd45haa5mh [ 2 ] ; real_T
lngqcivqqg [ 2 ] ; real_T ozmcvi5b31 [ 2 ] ; real_T hdadv1kifl [ 2 ] ; real_T
dhds3swxi0 [ 2 ] ; real_T oqc0b5does [ 2 ] ; real_T ieohtob00g [ 2 ] ; real_T
eiyeel4sy0 [ 2 ] ; real_T lgu0i43xpl [ 2 ] ; real_T n01cjpi1wr [ 2 ] ; real_T
bviufihvoa [ 2 ] ; real_T i0bfhtiky2 [ 2 ] ; real_T lq23ppo5t1 [ 2 ] ; real_T
ffetch5fn3 [ 2 ] ; real_T aj3ho3s0dk [ 2 ] ; real_T jjxvwteapb [ 2 ] ; real_T
mqctzxdg0o [ 2 ] ; real_T nnmjponkr5 [ 2 ] ; real_T ky3lis0a5h [ 2 ] ; real_T
iftqygn0yl [ 2 ] ; real_T ndljvxaegf [ 2 ] ; real_T k5hl2nqzw5 [ 2 ] ; real_T
exk4rrasy3 ; real_T oz0f0tgc3u ; void * ctycytv1sy ; void * m2wgpmvjvy ; void
* oqiguspxfn ; void * bigiskgqfj ; void * djxei35le3 ; void * dikp21y4fz ;
void * fqoysoyi51 ; void * dcyq5ws3qs ; void * ixn514xqve ; void * baw0vps5za
; struct { void * AQHandles ; } l4r0r31ubh ; void * bwezjkpzim ; void *
bi5anfwza5 ; void * gyisyfvd4m ; int_T h2bg1hotvh ; int_T glrxfbh4lh ; int8_T
oay4r2t35k ; int8_T itbalhm53p ; int8_T o3hz2jqijw ; int8_T a1ctm05tly ;
int8_T gib2w1lke4 ; int8_T ncvj32mtr4 ; int8_T ojwdvmmbzv ; int8_T blmbxswg4c
; int8_T eiqm0ndfyg ; int8_T pibyff0eql ; int8_T l2ec5jammp ; int8_T
iuadxosr40 ; int8_T iqaxcni1cl ; int8_T o4w412lyhz ; int8_T l4r1nqhb1g ;
int8_T gq3uhkidnv ; int8_T emuswyjmhs ; int8_T eqmbwqv4uu ; int8_T bxfixvxfdw
; int8_T dei1252yrl ; int8_T ogqlv3h0ke ; int8_T cpgojctm31 ; int8_T
hgiqtwygy1 ; int8_T ddu0m0coos ; int8_T plosltfnse ; int8_T i2wzmz255s ;
int8_T a1er4qtcuo ; int8_T chs3ptplyf ; int8_T jtymq0r3wm ; int8_T avqe5nr4ga
; int8_T h53r4zk0cj ; int8_T gktfawr2st ; int8_T bj4ejpvtiv ; int8_T
eyhr14ycpb ; int8_T i0q150xsvv ; int8_T inbjvmxdis ; int8_T lam1o4zpo4 ;
int8_T ot4rduhsnb ; int8_T gmlynrr0kl ; int8_T dfembtpulp ; int8_T mg55yfehmm
; int8_T knhfhd3fs2 ; boolean_T l4gngyaono ; boolean_T b0eijwagye ; boolean_T
oxke0tvvta ; boolean_T kr4a3cr0mf ; boolean_T gdblfyh2ib ; boolean_T
mmq5ttry5y ; boolean_T agwigalrld ; boolean_T fudyiowmus ; boolean_T
ftrjzsptcb ; boolean_T ciu1ngvzgz ; boolean_T csdsboi20g ; boolean_T
pxsktqyoas ; boolean_T bya2xiv5j3 ; boolean_T hnl0z3aet1 ; boolean_T
ecvugwr44v ; boolean_T dfdi1srria ; boolean_T ibwargsoog ; boolean_T
cb3tpkjx5q ; boolean_T fg4y55crvv ; boolean_T h4wxpouwz2 ; boolean_T
loncelqhha ; boolean_T afd0kmksnl ; boolean_T pv0jmcqzdg ; boolean_T
lbzvytb4d4 ; boolean_T pagreemrym ; boolean_T kkiofr3le1 ; boolean_T
pqnrnlrnbt ; boolean_T d1ahpk1qqr ; boolean_T hw5n3aggbh ; boolean_T
kgduzxmrs5 ; } DW ; typedef struct { real_T eouclruk0d [ 18 ] ; } X ; typedef
struct { real_T eouclruk0d [ 18 ] ; } XDot ; typedef struct { boolean_T
eouclruk0d [ 18 ] ; } XDis ; typedef struct { real_T eouclruk0d [ 18 ] ; }
CStateAbsTol ; typedef struct { real_T eouclruk0d [ 18 ] ; } CXPtMin ;
typedef struct { real_T eouclruk0d [ 18 ] ; } CXPtMax ; typedef struct {
real_T gxjau55nwb ; real_T n4sn2x5fn0 ; real_T nxyrsmdlot ; real_T jvkjfqj1ug
; real_T lge40yyb1l ; real_T bhqttm0cw1 ; real_T lphrnpyiuf ; real_T
js0gaa5lhi ; real_T eldmkhybny ; real_T j13ih0l24o ; real_T c0cdw0j0oq ;
real_T hii3txsdrd ; real_T ohmdhysjjv ; real_T bqtyhf5qky ; real_T iylzqyo3it
; real_T hxfby2z4m1 ; real_T gbmmjctguf ; real_T bkoti100ma ; real_T
ek2vlosxs0 ; real_T chg25xty1o ; real_T hrhezrd0m5 ; real_T gk4jezh4kc ;
real_T mccj2kywvs ; real_T njgojdbb1c ; real_T ipeaoktcgx ; real_T ijm0bfdiwn
; real_T osjkfrzf0x ; real_T myfa0e3bis ; } ZCV ; typedef struct {
rtwCAPI_ModelMappingInfo mmi ; } DataMapInfo ; struct P_ {
struct_C8RyKn0xmbw5jBmaos4hYE Fixed ; struct_0rpQxwwyqCy2sr7lBYKBtB SimParams
; real_T x [ 14 ] ; real_T NormalizeVector_maxzero ; real_T
NormalizeVector_maxzero_ab1bjfgz5c ; real_T
NormalizeVector_maxzero_jpxndxsk21 ; real_T NormalizeVector1_maxzero ; real_T
NormalizeVector_maxzero_hach1vuaw5 ; real_T
NormalizeVector_maxzero_nzh4nf0nvr ; real_T
NormalizeVector_maxzero_pjti1ch2qi ; real_T
NormalizeVector1_maxzero_hfyr12og40 ; real_T
NormalizeVector_maxzero_jvlgwkux2b ; real_T
NormalizeVector_maxzero_k3zfgywmre ; real_T
NormalizeVector_maxzero_jceokergvn ; real_T
NormalizeVector1_maxzero_edxcmuydlr ; real_T
NormalizeVector_maxzero_levzkkbbsj ; real_T
NormalizeVector_maxzero_kt51ehrppa ; real_T
NormalizeVector_maxzero_hbfzo5qsh2 ; real_T
NormalizeVector1_maxzero_os1qknjcxt ; real_T
NormalizeVector_maxzero_anlse3xi4h ; real_T
NormalizeVector_maxzero_dfin15pekf ; real_T
NormalizeVector_maxzero_nfiirapz00 ; real_T
NormalizeVector1_maxzero_pbowfqqdx5 ; real_T
NormalizeVector_maxzero_ap3h53mgxx ; real_T
NormalizeVector_maxzero_kytveykcqy ; real_T
NormalizeVector_maxzero_ihmtmp1vz1 ; real_T
NormalizeVector1_maxzero_byir5pwee4 ; real_T
NormalizeVector_maxzero_fqwaltgljo ; real_T
NormalizeVector_maxzero_hkbi1qbiut ; real_T
NormalizeVector_maxzero_dnlmszn2ri ; real_T
NormalizeVector1_maxzero_j5bp1hk5bd ; real_T Gain_Gain ; real_T
Gain_Gain_iby0aikld1 ; real_T Gain1_Gain ; real_T Gain_Gain_o3ussx2pey ;
real_T Gain1_Gain_inzuw1ftmc ; real_T Gain_Gain_ch2juginai ; real_T
Gain1_Gain_iijf52nd0p ; real_T Gain_Gain_cnz04octod ; real_T
Gain1_Gain_kgqmwgbxmo ; real_T Gain_Gain_mxfyehzrbo ; real_T
Gain1_Gain_pklggewqk2 ; real_T Gain_Gain_ch2qj2nra0 ; real_T
Gain1_Gain_mt5ivn2ky5 ; real_T Gain_Gain_kzgbyxe2u1 ; real_T
Gain1_Gain_fzh3lvxzcm ; real_T Gain_Gain_gnnygtgrhf ; real_T Constant_Value ;
real_T Constant_Value_foidw3bu40 ; real_T MagnetDipoleMoment_Value [ 3 ] ;
real_T Constant_Value_ffpsqzlbdy ; real_T Constant1_Value [ 9 ] ; real_T
Constant_Value_aldt0rz0ce ; real_T Constant_Value_hsgo25r4sk ; real_T
Constant_Value_lu1rr1djht ; real_T Constant_Value_hrhqyti5xa ; real_T
Constant_Value_bgbfcrykss ; real_T MagnetDipoleMoment_Value_hus3icj2zx [ 3 ]
; real_T Constant_Value_ccdbgnkv2m ; real_T Constant1_Value_mat3dcwilq [ 9 ]
; real_T Constant_Value_ag31g0q0wg ; real_T Constant_Value_lpo3xrcqjh ;
real_T Constant_Value_kebaex1duz ; real_T Constant_Value_crcwx2inuv ; real_T
Constant_Value_p4qexwk3wd ; real_T MagnetDipoleMoment_Value_m4h5ot0d2d [ 3 ]
; real_T Constant_Value_pxe2jyk15f ; real_T Constant1_Value_obh4dra0ft [ 9 ]
; real_T Constant_Value_a4fqp5ywdt ; real_T Constant_Value_lwhq5q2fck ;
real_T Constant_Value_mqmlcee2pu ; real_T Constant_Value_mulpe1a5eh ; real_T
Constant_Value_f4bgccfy2v ; real_T MagnetDipoleMoment_Value_am3fjqdljc [ 3 ]
; real_T Constant_Value_kyzvzubrmq ; real_T Constant1_Value_htaifltrf1 [ 9 ]
; real_T Constant_Value_hooldt5mbh ; real_T Constant_Value_fyrymlyvwp ;
real_T Constant_Value_iinhm22lbs ; real_T Constant_Value_o54qqmia1e ; real_T
Constant_Value_fkwr40m5ih ; real_T MagnetDipoleMoment_Value_hqiqmgkprz [ 3 ]
; real_T Constant_Value_c50rewkplm ; real_T Constant1_Value_ajigyq4ohn [ 9 ]
; real_T Constant_Value_gxy2qdgtgj ; real_T Constant_Value_oqtflizduh ;
real_T Constant_Value_e34xedeirm ; real_T Constant_Value_l5p2q1mni4 ; real_T
Constant_Value_cuvlenn5ze ; real_T MagnetDipoleMoment_Value_clffg1mnnq [ 3 ]
; real_T Constant_Value_pik4yifivl ; real_T Constant1_Value_dagdlneepy [ 9 ]
; real_T Constant_Value_ndmgh4xmik ; real_T Constant_Value_kqizw1abe3 ;
real_T Constant_Value_gimcsalf5s ; real_T Constant_Value_h2tkrf03re ; real_T
Constant_Value_oc3ljvgxhv ; real_T MagnetDipoleMoment_Value_af4j3ouf2p [ 3 ]
; real_T Constant_Value_nu3mnt2521 ; real_T Constant1_Value_k2adofcs3i [ 9 ]
; real_T Constant_Value_m0tx3vya3h ; real_T Constant_Value_jgubgdhktv ;
real_T Constant_Value_nfejvnig4t ; } ; extern const char_T *
RT_MEMORY_ALLOCATION_ERROR ; extern B rtB ; extern X rtX ; extern DW rtDW ;
extern P rtP ; extern mxArray * mr_MagneticBasket_Simscape_Optimizer_GetDWork
( ) ; extern void mr_MagneticBasket_Simscape_Optimizer_SetDWork ( const
mxArray * ssDW ) ; extern mxArray *
mr_MagneticBasket_Simscape_Optimizer_GetSimStateDisallowedBlocks ( ) ; extern
const rtwCAPI_ModelMappingStaticInfo *
MagneticBasket_Simscape_Optimizer_GetCAPIStaticMap ( void ) ; extern
SimStruct * const rtS ; extern DataMapInfo * rt_dataMapInfoPtr ; extern
rtwCAPI_ModelMappingInfo * rt_modelMapInfoPtr ; void MdlOutputs ( int_T tid )
; void MdlOutputsParameterSampleTime ( int_T tid ) ; void MdlUpdate ( int_T
tid ) ; void MdlTerminate ( void ) ; void MdlInitializeSizes ( void ) ; void
MdlInitializeSampleTimes ( void ) ; SimStruct * raccel_register_model (
ssExecutionInfo * executionInfo ) ;
#endif
