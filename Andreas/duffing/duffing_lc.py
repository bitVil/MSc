import auto


sysname  = "duffing_lc"
unames   = dict( enumerate( ["x1", "x2", "l1", "l2"], start=1) )
parnames = dict( enumerate( ["omega", "d", "f"], start=1) )

# system parameters in .c file
f = 0.0
d = 0.2
omega = 1.0

print(sysname + " is started")

# first increase f
f = 5
r1 = auto.run(
    e = sysname, 
    unames = unames,
    parnames = parnames,
    NDIM=   len(unames), IPS =  2, IRS =   0, ILP =   1,
    ICP =  ["f", "PERIOD"],
    NTST=  20, NCOL=   4, IAD =   3, ISP =   2, ISW = 1, IPLT= 9, NBC= 0, NINT= 0,
    NMX=  200, NPR=  100, MXBF=   0, IID =   2, ITMX= 8, ITNW= 7, NWTN= 3, JAC= 0,
    EPSL= 1e-07, EPSU = 1e-07, EPSS =0.0001,
    DS  =  0.05, DSMIN=  0.001, DSMAX=   0.1, IADS=   1,
    NPAR = len(parnames), THL =  {}, THU =  {},
    STOP = [],
    UZSTOP = {"f": f}
)


# auto.plot(r1("UZ"),solution_y = "x2")
# raw_input("Press Enter to continue...")

tmp = auto.run(r1("UZ1"),ICP=['omega'],UZSTOP = {"omega": 7}, NMX=3000)
r2 = auto.run(tmp("UZ1"), NTST=100, DS='-',UZSTOP = {"omega": 0.05})

r2 = r2.relabel()


sb1 = r2 + auto.run(r2("BP2"), ICP=['omega',"PERIOD"], ISW=-1,  NMX=200, NPR=10  )
sb1 = sb1.relabel()

auto.plot(sb1,stability=True,bifurcation_y="MAX x1")



tmp = auto.run(r2("LP2"),DS="-", ICP=['omega',"f","PERIOD"],UZSTOP = {"omega": 5, "f": 25}, NMX=3000, ISW=2, ILP=0)

r4 =  auto.run(tmp("EP"))
r4 = r4 + auto.run(tmp("EP"),DS="-")

tmp = auto.run(r2("BP1"),DS="-", ICP=['omega',"f","PERIOD"],UZSTOP = {"omega": 5, "f": 25}, NMX=3000, ISW=2, ILP=0)

r5 =  auto.run(tmp("EP"))
r5 =  r5 + auto.run(tmp("EP"),DS="-")

bifu = (r4 + r5).relabel()

namestr = "bifu_d"+ str(d)

fig = auto.plot(bifu, bifurcation_x = "omega", bifurcation_y="f",top_title=namestr,maxx=4)


fig.savefig(namestr+".pdf")

raw_input("Press Enter to continue...")
