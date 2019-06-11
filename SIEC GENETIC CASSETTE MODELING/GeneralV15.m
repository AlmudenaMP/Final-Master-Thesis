%AUTHORS: ALMUDENA MÉNDEZ PÉREZ & JESÚS DAZA GARCÍA

clear;
to = 0;
tf = 21600;

% Initial conditions

%v15 OFF STATE
%normalo = [0 0 0 0 0 0 0 0 0 0];

%v15 ON STATE
normalo = [0 133900 0 0 0 76000 0 0 0 0];

options = odeset('NonNegative',2);
[t, y] = ode45('v15',[to tf],normalo,options);
plot(t,y(:,2),'color','b'); hold on;
plot(t,y(:,4),'color','g'); hold on;
plot(t,y(:,6),'color','r'); hold on;
plot(t,y(:,8),'color','y'); hold on; 
plot(t,y(:,10),'color','m');

title('Protein concentrations')
xlabel('Time')
ylabel('Concentration')
legend('tetR','cI','LacI','TSS3 reporter','PROTEASE')

disp(y(2));
disp(y(4)); 
disp(y(6)); 
disp(y(8)); 
disp(y(10));
