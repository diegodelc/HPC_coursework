clear
close all

threads = [1,2,3,4,5,6,7,8,9,10];

runtimes = [3.626,2.156,1.577,1.281,1.027,0.871,0.787,0.864,0.916,1.040];

runtimesO3flag = [0.886,0.664,0.450,0.404,0.462,0.509,0.455,0.304,0.366,0.352];






%ideal = @(x) runtimes(1).*(0.5).^(x-1)

ideal = @(x) runtimes(1)./x;
idelO3flag = @(x) runtimesO3flag(1)./x;

figure()

loglog(threads,runtimes,'-ob')
hold on
loglog(ideal(threads),'-kx')
hold on
loglog(threads,runtimesO3flag,'-ob')
hold on
loglog(idelO3flag(threads),'-kx')
hold on
txt = '-O0 flag';
text(10^(0.1),10^(0.3),txt)
hold on
txt = '-O3 flag';
text(10^(0.1),10^(-0.3),txt)
legend("Scaling found","Ideal scaling")
xlabel("log(thread count)")
ylabel("log(runtime)")