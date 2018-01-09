%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sahebeh Dadboud: 1569395
%Assignment3-exe3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% advect - Program to solve the advection equation 
% using the various hyperbolic PDE schemes
clear all; % help advect;  % Clear memory and print header

%% * Select numerical parameters (time step, grid spacing, etc.).
method = menu('Choose a numerical method:', ...
       'FTCS','Lax','Lax-Wendroff');
N = input('Enter number of grid points: ');
L = 1.;     % System size
h = L/(N-1);    % Grid spacing for Dirichlet BC
c = 1;      % Wave speed
fprintf('Time for wave to move one grid spacing is %g\n',h/c);
tau = input('Enter time step: ');
coeff = -c*tau/(2.*h);  % Coefficient used by all schemes
coefflw = 2*coeff^2;    % Coefficient used by L-W scheme
fprintf('Wave circles system in %g steps\n',L/(c*tau));
nStep = input('Enter number of steps: ');

%% * Set initial and boundary conditions (Dirichlet)

x = (1:N)*h - L/2 - h;  % Coordinates of grid points
w = 10*pi;
% set scalar field a  to zero.
a(1:N) = 0;


%% * Initialize plotting variables.
iplot = 1;          % Plot counter
aplot(:,1) = a(:);  % Record the initial state
tplot(1) = 0;       % Record the initial time (t=0)
nplots = 50;        % Desired number of plots
plotStep = nStep/nplots; % Number of steps between plots

%% * Loop over desired number of steps.
for iStep=1:nStep  %% MAIN LOOP %%
  t = tau*iStep;
  a_new(1:N) = 0;
  % Dirichlet boundary conditions:
  a_new(1) = sin(w*t);
  a_new(N) = 0; 
  %* Compute new values of wave amplitude using FTCS, 
  %  Lax or Lax-Wendroff method.
  if( method == 1 )      %%% FTCS method %%% a(1:N) = a(1:N) + coeff*(a(ip)-a(im)); 
    for i = 2:(N-1)
      a_new(i) = a(i) + coeff*(a(i+1)-a(i-1));
    end
  elseif( method == 2 )  %%% Lax method %%% a(1:N) = .5*(a(ip)+a(im)) + coeff*(a(ip)-a(im)); 
    for i=2:(N-1)
      a_new(i) = .5*(a(i+1)+a(i-1)) + coeff*(a(i+1)-a(i-1));
    end
  else  %%% Lax-Wendroff method %%% a(1:N) = a(1:N) + coeff*(a(ip)-a(im)) + coefflw*(a(ip)+a(im)-2*a(1:N)); 
    for i=2:(N-1)
      a_new(i) = a(i) + coeff*(a(i+1)-a(i-1)) + coefflw*(a(i+1)+a(i-1)-2*a(i)); 
    end
  end   
  a(:) = a_new(:);

  %* Periodically record a(t) for plotting.
  if( rem(iStep,plotStep) < 1 )  % Every plot_iter steps record 
    iplot = iplot+1;
    aplot(:,iplot) = a(:);       % Record a(i) for ploting
    tplot(iplot) = tau*iStep;
    fprintf('%g out of %g steps completed\n',iStep,nStep);
  end
end

%% * Plot the initial and final states.
figure(1); clf;  % Clear figure 1 window and bring forward
plot(x,aplot(:,1),'-',x,a,'--');
legend('Initial  ','Final');
xlabel('x');  ylabel('a(x,t)');
pause(1);    % Pause 1 second between plots

%% * Plot the wave amplitude versus position and time
figure(2); clf;  % Clear figure 2 window and bring forward
mesh(tplot,x,aplot);
ylabel('Position');  xlabel('Time'); zlabel('Amplitude');
view([-70 50]);  % Better view from this angle