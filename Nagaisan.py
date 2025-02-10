import numpy as np
from scipy import interpolate
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

def main():
    with open("/Users/shohtatakami/physics/COMETDS/Codes/Nagaisan/Y=100_FieldMap.csv") as f:
        data_diff = []
        for line in f.readlines():
            print(line)
            for i in line.split("\n")[0].split("\t"):
                data_diff.append(float(i))
                #data_diff.append([float(i) for i in line.split("\n")[1].split("\t")])
        

    with open("/Users/shohtatakami/physics/COMETDS/Codes/Nagaisan/Y=0_FieldMap.csv") as f:
        data = []
        for line in f.readlines():
            data.append([float(i) for i in line.split("\n")[0].split("\t")])

    pdf = PdfPages("using3DMapFile.pdf")
    xz_data = []
    zlist = []
    for i in range(len(data)):
        if i==0:
            x = data[i][0] 
            linez = []
        elif x != data[i][0]:
            xz_data.append(linez)
            x = data[i][0]
            linez = []
            zlist = []
        linez.append(data[i][6])
        zlist.append(data[i][2])

    plt.imshow(xz_data)
    pdf.savefig()
    plt.clf()

    z = []
    x0_data = []
    for i in range(len(data)):
        if data[i][0] != 0:
            continue
        z.append(data[i][2])
        x0_data.append(data[i][6])

    diff = 10 #mm
    zlim = []
    #xlist = np.linspace(-800,800,17)
    xlist = [-760,-500,0,500,760]
    #xlist = [760,500,0,-500,-760]
    corrB = [[] for i in range(len(xlist))]
    origB = [[] for i in range(len(xlist))]
    for z in zlist:
        if np.abs(z) > 1000:
            continue
        zlim.append(z)
        tantheta = diff/1500
        zsorted = []
        xdata = []
        for element in data:
            if element[2]==z:
                zsorted.append(element[6])
                xdata.append(element[0])
        #print(zsorted, xdata)
        f = interpolate.interp1d(xdata, zsorted,kind='linear')
        fullx = np.linspace(-1000,1000,1001)
        fullB = f(fullx)
        corrx = z * tantheta
        print(corrx)
        for i in range(len(xlist)):
            corrB[i].append(f(xlist[i]+corrx))
            origB[i].append(f(xlist[i]))
        print(np.array(corrB[i])- np.array(origB[i]))


    def pol2(x,a,b,c):
        return a*(x-b)**2 + c

    peakm_orig = []
    peakm_corr = []
    peake_orig = []
    peake_corr = []
    for i in range(len(xlist)):
        z_prec = np.linspace(-1000,1000,1001)
        fit_range = np.array([abs(z) for z in zlim])<250
        plt.plot(zlim,origB[i],'o')
        plt.plot(zlim,corrB[i],'o')
        popt_orig,pcov_orig = curve_fit(pol2, np.array(zlim)[fit_range], np.array(origB[i])[fit_range])
        popt_corr,pcov_corr = curve_fit(pol2, np.array(zlim)[fit_range], np.array(corrB[i])[fit_range])
        perr_orig = np.sqrt(np.diag(pcov_orig))
        perr_corr = np.sqrt(np.diag(pcov_corr))
        plt.plot(z_prec,pol2(z_prec,*popt_orig), color='tab:blue')
        plt.plot(z_prec,pol2(z_prec,*popt_corr), color='tab:orange')
        plt.axvline(popt_orig[1],linestyle=':')
        plt.axvline(popt_corr[1],linestyle=':',color='tab:orange')
        plt.title(f"{xlist[i]:.2f}")
        filename = f"x_{xlist[i]:.0f}_diff_{diff:.0f}.csv"  # ファイル名に xlist[i] と diff を含める
        df = pd.DataFrame({"z": zlim, "Corrected B": corrB[i], "Original B": origB[i]})
        df.to_csv(filename, index=False)  # index=False で余計なインデックスを出力しない
        print(f"Saved: {filename}")
        pdf.savefig()
        plt.show(block=False)
        input("Enter を押すと次へ進みます...")  # ユーザーの入力待ち
        plt.clf()
        peakm_orig.append(popt_orig[1])
        peakm_corr.append(popt_corr[1])
        peake_orig.append(perr_orig[1])
        peake_corr.append(perr_corr[1])
    
    plt.errorbar(xlist, peakm_orig, yerr=peake_orig, fmt='o',label='No inclination')
    plt.errorbar(xlist, peakm_corr, yerr=peake_corr, fmt='o',label='w/ inclination')
    plt.axvline(0,color='gray',linestyle=':')
    plt.axhline(0,color='gray',linestyle=':')
    plt.title(f'diff = {diff} mm at the edge')
    plt.xlabel('x [mm]')
    plt.ylabel('peak position of z [mm]')
    plt.legend()
    pdf.savefig()
    plt.clf()

    pdf.close()

if __name__=='__main__':
    main()
