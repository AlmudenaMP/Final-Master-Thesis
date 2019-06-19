%AUTHORS: ALMUDENA MÉNDEZ PÉREZ & JESÚS DAZA GARCÍA

clear;
to = 0;
tf = 21600;

%OFF STATE 

%Initial conditions for v1 v3 v10
%normalo = [0 0 0 0 0 0 0 0];

%ON STATE

%Initial conditions for v1 
%normalo = [0 133900 0 0 0 75900 0 0];

%Initial conditions for v3 
%normalo = [0 133900 0 0 0 0 0 86780];

%Initial conditions for v10 
%normalo = [0 133900 0 0 0 759.8 0 0];
 
options = odeset('NonNegative',2);
[t, y] = ode45('v10',[to tf],normalo,options);
plot(t,y(:,2),'color','b'); hold on;
plot(t,y(:,4),'color','g'); hold on;
plot(t,y(:,6),'color','r'); hold on;
plot(t,y(:,8),'color','y'); 

%title('Protein concentrations')
xlabel('Time (h.)')
ylabel('Protein concentrations')
legend('tetR','cI','LacI','T3SS reporter')

disp(y(2));
disp(y(4)); 
disp(y(6)); 
disp(y(8)); 