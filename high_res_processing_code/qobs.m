
 % output qobs profile from Neale and Hoskins to a file

 clear
for res = [0,4,8,16,32]

if res == 4
 dy = 48e3; % m
 ny = 360;
elseif res == 8
 dy = 96e3; % m
 ny = 180;
elseif res == 16
 dy = 192e3; % m
 ny = 90;

elseif res ==32 
 dy = 384e3; % m
 ny = 45;

elseif res ==0
 dy = 384e3; % m
 ny = 48;
else
sprintf('chose a wrong coarse factor');
end

% dy = 240e3; % m
% ny = 80;

 for j=1:ny
  lat(j)=dy*(j-0.5-ny/2)*360/4e7; % consistent with SAM expression
 end
 lat(1)
 lat(ny)

 polar = find(abs(lat)*pi/180>pi/3);
 
 zero_C = 273.15;

 control_sst = 27.0*(1-sind(3.0*lat/2.0).^2) + zero_C;
 control_sst(polar)=zero_C;

 flat_sst = 27.0*(1-sind(3.0*lat/2.0).^4) + zero_C;
 flat_sst(polar)=zero_C;

 qobs_sst = 0.5*(control_sst+flat_sst);

 plot(lat, flat_sst, lat, control_sst, lat, qobs_sst)

 filename = sprintf('ssty_%i',res)
 fileID = fopen(filename,'w'); 
 %fileID = fopen('ssty','w');
 for j=1:ny
  fprintf(fileID,'%6.2f\n',qobs_sst(j));
 end
 fclose(fileID);
 filename4k = sprintf('ssty4K_%i',res)
 fileID = fopen(filename4k,'w');
 %fileID = fopen('ssty_4K','w');
 for j=1:ny
  fprintf(fileID,'%6.2f\n',qobs_sst(j)+4);
 end
 fclose(fileID);
end
