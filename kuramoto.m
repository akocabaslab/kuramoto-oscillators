%2D Kuamoto

nx = 100;
ny = 100;

theta=rand(nx,ny) * 2 * pi - pi;


nt = 1000000;
omega = ones(nx,ny);



range = 1;


K0 =0.25;
b = 0.3*1.0000*ones(nx,ny);
alpha = 0.2*pi*ones(nx,ny);


% 
% 
mask =ones(nx, ny);




dt = 0.001;


for t = 1:1:nt
    % Loop over oscillators and update their phases
   
    neighbors = zeros(nx,ny);
    count=0;
 
     

       i=1;
       j=1;        
               
    right_neighbor=circshift(theta,[i 0]);
    left_neighbor=circshift(theta,[-i 0]);
    up_neighbor=circshift(theta,[0 j]);
    down_neighbor=circshift(theta,[0 -j]);
   
    %No-flux Boundary
   neighbors = K0*sin(right_neighbor-theta-alpha) + K0*sin(left_neighbor-theta-alpha) + K0*sin(up_neighbor-theta-alpha) + K0*sin(down_neighbor-theta-alpha) ;
 
    
   neighbors(1,2:ny-1)=2*K0*sin(left_neighbor(1,2:ny-1)-theta(1,2:ny-1)-alpha(1,2:ny-1))+K0*sin(up_neighbor(1,2:ny-1)-theta(1,2:ny-1)-alpha(1,2:ny-1)) + K0*sin(down_neighbor(1,2:ny-1)-theta(1,2:ny-1)-alpha(1,2:ny-1)) ;
    neighbors(2:nx-1,1)=K0*sin(right_neighbor(2:nx-1,1)-theta(2:nx-1,1)-alpha(2:nx-1,1)) + K0*sin(left_neighbor(2:nx-1,1)-theta(2:nx-1,1)-alpha(2:nx-1,1))+ 2*K0*sin(down_neighbor(2:nx-1,1)-theta(2:nx-1,1)-alpha(2:nx-1,1));
   neighbors(nx,2:ny-1)=2*K0*sin(right_neighbor(nx,2:ny-1)-theta(nx,2:ny-1)-alpha(nx,2:ny-1))+K0*sin(up_neighbor(nx,2:ny-1)-theta(nx,2:ny-1)-alpha(nx,2:ny-1)) + K0*sin(down_neighbor(nx,2:ny-1)-theta(nx,2:ny-1)-alpha(nx,2:ny-1)); 
   neighbors(2:nx-1,ny)=K0*sin(right_neighbor(2:nx-1,ny)-theta(2:nx-1,ny)-alpha(2:nx-1,ny)) + K0*sin(left_neighbor(2:nx-1,ny)-theta(2:nx-1,ny)-alpha(2:nx-1,ny))+ 2*K0*sin(up_neighbor(2:nx-1,ny)-theta(2:nx-1,ny)-alpha(2:nx-1,ny)) ;
    
    neighbors(1,1)=2*K0*sin(left_neighbor(1,1)-theta(1,1)-alpha(1,1))+2*K0*sin(down_neighbor(1,1)-theta(1,1)-alpha(1,1));
   neighbors(1,ny)=2*K0*sin(right_neighbor(1,ny)-theta(1,ny)-alpha(1,ny))+2*K0*sin(up_neighbor(1,ny)-theta(1,ny)-alpha(1,ny));
    neighbors(nx,1)=2*K0*sin(right_neighbor(nx,1)-theta(nx,1)-alpha(nx,1))+2*K0*sin(down_neighbor(nx,1)-theta(nx,1)-alpha(nx,1));
    neighbors(nx,ny)=2*K0*sin(left_neighbor(nx,ny)-theta(nx,ny)-alpha(nx,ny))+2*K0*sin(up_neighbor(nx,ny)-theta(nx,ny)-alpha(nx,ny));
   
   
   

   
   
   theta=theta+(1*omega-mask.*b.*sin(theta))*dt+mask.*neighbors*dt;
   theta=wrapToPi(theta);
   



    if mod(t,500) == 0
        

    figure(1);
    
    %Create a large canvas for combined images
    set(gcf, 'Position', [50, 100, 900, 900]); % Set figure size
    set(gcf, 'Visible', 'on'); % Ensure the figure is visible
    
    subplot(1,1, 1);
    imagesc(imresize((1 + sin(theta)) * 0.5, 3));
    colormap(gca, 'jet'); % Apply jet colormap only to this subplot
    colorbar; % Separate colorbar for this subplot
    title('Theta');
    axis image;


    end


end


