 % Plotting Vector Field, i.e. principle stresses
 close all
 clear all
 plot_direction_data;
 % Make plot
 scale_factor = 0.4;

 % Define characteristic length
 clx = 0;   
 cly = 0;
 for e = 1:size(IX,1)
     for i = 1:4
         for j = i+1:4
             clx = max(abs(  X(IX(e,i),1) - X(IX(e,j),1)   ),clx);
             cly = max(abs(  X(IX(e,i),2) - X(IX(e,j),2)   ),cly);
         end
     end
 end
 clmax = max(clx, cly);
 scal = clmax * scale_factor;
 % Make plot
 figure
 hold on
 for e = 1:size(IX,1)
     xx = X(IX(e,1:4),1);
     yy = X(IX(e,1:4),2);
     patch(xx,yy,[1 1 1]);
     % Find approx. center for isoparametric element
     xc = sum(X(IX(e,:),1))/size(IX,2);
     yc = sum(X(IX(e,:),2))/size(IX,2);
     % Directions for vect(:)
     vec = [cos(-vect(e)) sin(-vect(e)) ...
     cos(-vect(e)+pi/2) sin(-vect(e)+pi/2)];

     % Plot magnitude and direction
     cc = 'k';
     quiver(xc,yc,vec(1),vec(2),  scal,cc)
     quiver(xc,yc,-vec(1),-vec(2),scal,cc)
 end   
 title( 'Direction')
 axis equal;  axis off;  hold off
 set(gcf,'color',[ 1  1  1]);
