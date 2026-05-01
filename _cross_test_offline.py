import math, random, csv

THETA0=58.65; TOL_ANG=4.0
POLES_DEF=[
    (305., 25.,'D'),(125.,-25.,'A'),(35., 25.,'S1'),
    (215.,-25.,'S2'),(215., 25.,'S3'),(35.,-25.,'S4'),
]

def lb_xyz(l,b):
    lr,br=math.radians(l),math.radians(b)
    return (math.cos(br)*math.cos(lr), math.cos(br)*math.sin(lr), math.sin(br))

def radec_to_lb(ra,dec):
    ra_r=math.radians(ra); dc_r=math.radians(dec)
    RN=math.radians(192.85948); DN=math.radians(27.12825)
    b=math.asin(max(-1.,min(1.,
       math.sin(dc_r)*math.sin(DN)+math.cos(dc_r)*math.cos(DN)*math.cos(ra_r-RN))))
    y=math.cos(dc_r)*math.sin(ra_r-RN)
    x=math.cos(dc_r)*math.sin(DN)*math.cos(ra_r-RN)-math.sin(dc_r)*math.cos(DN)
    l=(math.degrees(math.atan2(y,x))+122.93192)%360.
    return l,math.degrees(b)

def ang_sep_xyz(v1,v2):
    return math.degrees(math.acos(max(-1.,min(1.,sum(a*b for a,b in zip(v1,v2))))))

def crossing_density(l,b):
    sv=lb_xyz(l,b); nc,which=0,[]
    for pl,pb,name in POLES_DEF:
        pv=lb_xyz(pl,pb)
        th=ang_sep_xyz(pv,sv)
        for n in [1,2,3]:
            if abs(th-n*THETA0)<TOL_ANG:
                nc+=1; which.append((name,n,th))
    return nc,which

def crossproduct_pa(ra_src,dec_src,which):
    """Best-pair cross product projected onto sky."""
    seen,pvs=[],[]
    for pname,n,theta in which:
        if pname in [x[1] for x in pvs]:
            continue
        for pl_l,pl_b,p in POLES_DEF:
            if p==pname:
                pvs.append((lb_xyz(pl_l,pl_b),pname)); break
    if len(pvs)<2:
        return None,None
    best_mag=0.; bcx=bcy=bcz=0.; bn1=bn2=''
    for i in range(len(pvs)):
        for j in range(i+1,len(pvs)):
            v1,n1=pvs[i]; v2,n2=pvs[j]
            cx=v1[1]*v2[2]-v1[2]*v2[1]
            cy=v1[2]*v2[0]-v1[0]*v2[2]
            cz=v1[0]*v2[1]-v1[1]*v2[0]
            mag=math.sqrt(cx*cx+cy*cy+cz*cz)
            if mag>best_mag:
                best_mag=mag; bcx,bcy,bcz=cx,cy,cz; bn1,bn2=n1,n2
    if best_mag<0.1:
        return None,None  # alle Paare fast antipodal
    bcx/=best_mag; bcy/=best_mag; bcz/=best_mag
    ra_r=math.radians(ra_src); dec_r=math.radians(dec_src)
    nh=(-math.sin(dec_r)*math.cos(ra_r),-math.sin(dec_r)*math.sin(ra_r),math.cos(dec_r))
    eh=(-math.sin(ra_r),math.cos(ra_r),0.)
    pn=bcx*nh[0]+bcy*nh[1]+bcz*nh[2]
    pe=bcx*eh[0]+bcy*eh[1]+bcz*eh[2]
    if abs(pn)<1e-10 and abs(pe)<1e-10:
        return None,None
    return math.degrees(math.atan2(pe,pn))%180, f'{bn1}x{bn2}'

def align(j_pa,e_pa):
    d=abs((j_pa%180)-(e_pa%180))
    return 180.-d if d>90. else d

# Literatur Jet-PA (mod 180 = Achse)
Lit={
    'M 87':291,'NGC 4486':291,'M 84':308,'NGC 4374':308,'NGC 4261':88,
    'NGC 4278':116,'NGC 4552':110,'M 89':110,'NGC 4649':83,'NGC 4636':142,
    'NGC 1275':160,'NGC 1052':65,'NGC 5128':51,'NGC 6251':251,'NGC 315':175,
    'NGC 1600':68,'NGC 3379':155,'NGC 4459':90,'NGC 4473':92,'NGC 4564':47,
    'NGC 4697':65,'NGC 4742':5,'NGC 4889':78,'NGC 7052':100,'NGC 821':50,
    'NGC 383':15,'NGC 1265':125,'NGC 3115':50,'NGC 3608':100,'NGC 4168':36,
    'NGC 7768':30,'NGC 5845':81,
}

rows=[]
with open(r'c:\Users\kalle\Neuer Ordner (4)\results\catalogs\smbh_extended.csv',
          newline='', encoding='utf-8') as f:
    for r in csv.DictReader(f):
        try: ra=float(r['RA_deg']); dec=float(r['Dec_deg'])
        except: continue
        l,b=radec_to_lb(ra,dec); nc,which=crossing_density(l,b)
        rows.append({'name':r['Name'].strip(),'ra':ra,'dec':dec,'nc':nc,'which':which})

jets=[]
for s in rows:
    nm=s['name']
    for k,v in Lit.items():
        if k.upper().replace(' ','')==nm.upper().replace(' ',''):
            cp,cpn=crossproduct_pa(s['ra'],s['dec'],s['which'])
            jets.append({'name':nm,'nc':s['nc'],'which':s['which'],
                         'jet_pa':v%180,'cross_pa':cp,'cross_poles':cpn,
                         'cross_delta':align(v,cp) if cp is not None else None})
            break

print('=== Kreuzprodukt-Test: Fragmentierungs-Modell (Lit-Quellen) ===')
print('Prediction: Jet-Achse || P1 x P2 (freie Richtung bei Schalen-Kollision)')
print()

nc2=[j for j in jets if j['nc']>=2 and j['cross_delta'] is not None]
nc0=[j for j in jets if j['nc']==0 and j['cross_delta'] is not None]
nc1=[j for j in jets if j['nc']==1 and j['cross_delta'] is not None]
all_v=[j for j in jets if j['cross_delta'] is not None]
nc2_antipodal=[j for j in jets if j['nc']>=2 and j['cross_pa'] is None]

def stats(lst,lbl):
    if not lst: print(f'  {lbl}: N=0'); return
    n=len(lst); ds=sorted(j['cross_delta'] for j in lst)
    nal=sum(1 for d in ds if d<30)
    print(f'  {lbl:15s} N={n:3d}  aligned<30deg={nal}/{n}({100*nal/n:.0f}%)'
          f'  mean={sum(ds)/n:.1f}  med={ds[n//2]:.1f}')

stats(nc0,'nc=0 (ref)')
stats(nc1,'nc=1 (ref)')
stats(nc2,'nc>=2 (test)')
stats(all_v,'alle')
print(f'  Zufallserwartung: 33.3%')
print(f'  nc>=2 ohne Prediction (antipodal Paare): N={len(nc2_antipodal)}')
if nc2_antipodal:
    for j in nc2_antipodal:
        rings=','.join(f'{p}n{n}' for p,n,_ in j['which'])
        print(f'    {j["name"]}  nc={j["nc"]}  ringe={rings}')

print()
# MC
random.seed(42)
N=50000
mc_h=sum(1 for _ in range(N) if align(random.uniform(0,180),random.uniform(0,180))<30)
fmc=mc_h/N
print(f'MC-Referenz: {100*fmc:.1f}% (N={N})')
if nc2:
    fobs=sum(1 for j in nc2 if j['cross_delta']<30)/len(nc2)
    ratio=fobs/fmc if fmc>0 else float('nan')
    flag='*** SIGNAL ***' if ratio>1.5 else ('kein Signal' if ratio<0.8 else 'schwach')
    print(f'nc>=2 Kreuzprodukt: obs={100*fobs:.1f}%  MC={100*fmc:.1f}%  Ratio={ratio:.2f}  {flag}')

print()
print('nc>=2 Quellliste (nach cross_delta sortiert):')
print(f'  {"Name":20s} nc  {"PaarXP":8s}  JetPA  P1xP2  delta  (Ringe)')
for j in sorted([x for x in jets if x['nc']>=2], key=lambda x:(x['cross_delta'] is None, x['cross_delta'] or 999)):
    rings=','.join(f'{p}n{n}' for p,n,_ in j['which'])
    if j['cross_pa'] is None:
        print(f'  {j["name"]:20s} {j["nc"]:2d}  {"antipodal":8s}  {j["jet_pa"]:5.1f}  {"---":>6}  {"---":>5}  ({rings})')
    else:
        mk='ok' if j['cross_delta']<30 else '  '
        print(f'  {j["name"]:20s} {j["nc"]:2d}  {j["cross_poles"]:8s}  {j["jet_pa"]:5.1f}  {j["cross_pa"]:5.1f}  {j["cross_delta"]:5.1f} {mk}  ({rings})')
