function mkmaskGBT(galname,version,dilation)

% Adapted from code originally written by Alberto Bolatto.

% Uses galaxy parameter file galaxy_mask.csv
% Methods:
% - Havfield : uses the Ha intensity and velocity field in CALIFA and
%              Vicente's relation between linewidth and distance to the
%              center to produce mask.
% - rotcur   : uses the Ha intensity and a standard flat rotation curve
%              at vmaxg km/s and produces two masks for 180 deg rotation
%              orientations with the geometry provided. Intended for cases
%              where the Ha velocities are bad.
% - flat     : masks using an Ha intensity threshold and -/+ 2*vmaxg km/s
% - block    : masks the R25 region over -/+ 2*vmaxg km/s


snrvcutoff=500; % funky V SNR in CALIFA cubes
snrsfrcutoff=6;
vchw=10; % in km/s
vspan=1000; % range of velocities for mask cube

gpar=parseCsvTable('galaxy_mask.csv');
% This was designed to do one galaxy at a time, but I twisted the logic to
% enable doing all galaxies in the mask list if no parameters are passed.
% It's a bit contrived, but it works.
ngal=length(gpar.galname);
if nargin
    ngal=1;
end
if (nargin<3)
    dilation=2;
end
se=strel('cuboid',[dilation*2+1,dilation*2+1,dilation+1]);
nn=0;
while nn<ngal
    nn=nn+1;
    if (ngal>1)
        galname=gpar.galname{nn};
        version=gpar.version(nn);
    end
    fprintf(1,'%s, version %d\n',galname,version);
    ix=find(strcmp(galname,gpar.galname));
    if isempty(ix)
        fprintf(1,'Galaxy not found. Need parameters in galaxy_mask.csv\n');
        return;
    else
        if (length(ix)>1)
            ixx=find(gpar.version(ix)==version);
            ix=ix(ixx);
        else
            version=1;
        end
        logd25=gpar.logd25(ix);
        vmaxg=gpar.vmaxg(ix);
        PA=gpar.PA(ix);
        INC=gpar.INC(ix);
        maskat=gpar.maskat(ix);
        vmethod=gpar.method{ix};
    end
    
    R25=(0.1*60/2)*10^logd25;
    PA=90-PA;
    R=[0,R25/10,80];
    velR=[0,vmaxg,vmaxg];
    ddir='./';
    fname=[galname '.Pipe3D.cube.fits'];
    if ~exist(fname,'file')
        fprintf(1,'Pipe3D file not present. Downloading and decompressing...\n');
        websave([fname '.gz'],['https://ifs.astroscu.unam.mx/CALIFA/V500/v2.3/pyPipe3D/' fname '.gz']);
        system(['gunzip ' fname '.gz']);
        fprintf(1,'Done!\n');
    end
    
    c=2.99792458e5; % speed of light
    p=parseGBTCatalog('GBTEDGE.cat');
    RA=p.RA(strcmp(galname,p.NAME));
    DEC=p.DEC(strcmp(galname,p.NAME));
    vgal=c*str2double(string(p.Z(strcmp(galname,p.NAME))));
    
    %r=rfits([ddir fname],'xten3'); % OLD pipe3D setting (for NGC2596)
    r=rfits([ddir fname],'xten5');
    u=permute(r.data,[2 1 3]);

    % OLD pipe3D indices (for NGC2596)
    %Hae=u(:,:,250);
    %Hav=u(:,:,97);
    %Have=u(:,:,301);

    % NEW pipe3D indices
    Hae=u(:,:,262);    
    Hav=u(:,:,100);    
    Have=u(:,:,316);

    Ha=u(:,:,46); 
    Hb=u(:,:,29);
    Ha(Ha<0)=0;
    Hb(Hb<0)=0;
    
    % computation of extinction and CO intensity
    Aha=5.86*log10(Ha./(2.86*Hb)); % extinction correction computation
    Aha(Aha>3)=3;
    Aha(Aha<0)=0;
    DistMpc=1; % Dummy distance, just simplifies the calculations
    sfr=1e-16*5.5e-42*Ha.*10.^(Aha/2.5)*4*pi*(DistMpc*1e6*3.09e18)^2;
    snr=Ha./Hae;
    snrv=Hav./Have;
    sfr(isnan(sfr))=0;
    sfr(snr<snrsfrcutoff)=0; % cleaning up a bit by imposing some SNR cutoff
    Mmol=sfr*1.3e9; % computing Mmol from SFR with 1.3 Gyr timescale
    ScoDv=Mmol/(4.36*1.05e4*DistMpc^2); % computing Sco*dv in Jy km/s equivalent to CO mass
    gain=3.2e-4*(100)^2; % 100m dish gain in K/Jy
    comap=ScoDv*gain; % converting to K km/s (in Tmb)
    
    x=r.x{1}-RA*15;
    ra=x(1,:)*3600;
    y=r.x{2}-DEC;
    de=y(:,1)*3600;
    comapconv=convolve(x*3600,y*3600,comap,2.7,6,1);
    [a,b]=size(x);
    rot=[cosd(PA) sind(PA);-sind(PA) cosd(PA)];
    xs=reshape(x,[1,a*b]);
    ys=reshape(y,[1,a*b]);
    xy=[xs;ys];
    xyr=rot*xy;
    xr=reshape(xyr(1,:),[a,b]);
    yr=reshape(xyr(2,:)/cosd(INC),[a,b]);
    rr25=sqrt((xr*3600).^2+(yr*3600).^2)/R25;
    fwhm=[0,300;0.15,200;0.4,100;0.6,70;0.8,70;1,70;1.2,70;10,70];
    vwid=spline(fwhm(:,1),fwhm(:,2),rr25);
    vcen=Hav-vgal; vcen(snrv<snrvcutoff)=NaN; vcen(abs(vcen)>500)=NaN;
    nchan=1+round(vspan/vchw);
    vel=(-vspan/2):vchw:(vspan/2);
    z=(vel+vgal)/c; velo=c*z./(1+z);
    %vop=velo./(1-velo/c); vel=vop-vgal;
    %%cube=zeros(b,a,nchan);
    [~,~,vv]=meshgrid(ra,de,vel);
    
    clear im
    im.naxis=3;
    im.bitpix=-32;
    im.numpt=r.numpt;
    im.numpt(3)=length(vel);
    im.crval=r.crval;
    im.crval(3)=(velo(fix(nchan/2)+1));
    im.crpix=r.crpix;
    im.crpix(3)=fix(nchan/2)+1;
    im.ctype=r.ctype;
    im.ctype{3}='VRAD';
    im.cdelt(1)=r.cd1_1;
    im.cdelt(2)=r.cd2_2;
    im.cdelt(3)=-mean(diff(velo)); % Flipping 3rd axis order
    %im.cd1_1=r.cd1_1;
    %im.cd2_2=r.cd2_2;
    %im.cd2_1=0;
    %im.cd1_2=0;
    %im.cd3_3=im.cdelt(3);
    im.cunit=r.cunit;
    im.cunit{3}='km s-1';
    im.equinox=2000;
    im.coormode='ABSCOOR';
    im.restfreq=115271.202e6;
    im.specsys='LSRK';
    im.ssysobs='TOPOCENT';
    im.radesys='FK5';
    im.bmaj=2.7/3600;
    im.bmin=2.7/3600;
    im.bpa=0;
    
    if strcmp(vmethod,'rotcur') || strcmp(vmethod,'rotnoHa')
        vcen=mkvelfield(r,(r.crval(1)-RA*15)*cosd(r.crval(2))*3600,(r.crval(2)-DEC)*3600,PA,INC,R,velR);
    end

    vfac=2; %factor to increase vmaxg width
    
    % dv=(vv-vcen);
    % cube=comapconv.*exp(-dv.^2./(vwid/2.35).^2);
    % cubemask=zeros(size(cube));
    % cubemask(cube>maskat)=1;
    % im.data=permute(cubemask,[2 1 3]);
    
    switch vmethod
        case 'rotcur'
            dv=(vv-vcen);
            cube=comapconv.*exp(-dv.^2./(2*(vwid/2.35).^2));
            cubemask=zeros(size(cube));
            cubemask(cube>maskat)=1;
            %im.data=permute(cubemask,[2 1 3]);
            im.data=flip(permute(cubemask,[2 1 3]), 3);
            wfits(im,['masks/from_matlab/mask_' galname '_' vmethod '+_v' int2str(version) '.fits']);
            vcen=mkvelfield(r,(r.crval(1)-RA*15)*cosd(r.crval(2))*3600,(r.crval(2)-DEC)*3600,PA+180,INC,R,velR);
            dv=(vv-vcen);
            cube=comapconv.*exp(-dv.^2./(2*(vwid/2.35).^2));
            cubemask=zeros(size(cube));
            cubemask(cube>maskat)=1;
            im.data=flip(permute(cubemask,[2 1 3]),3);
            %im.data=permute(cubemask,[2 1 3]);
            if dilation>0, im.data=imdilate(im.data,se); end
            wfits(im,['masks/from_matlab/mask_' galname '_' vmethod '-_v' int2str(version) '.fits']);
        case 'Havfield'
            dv=(vv-vcen);
            cube=comapconv.*exp(-dv.^2./(2*(vwid/2.35).^2));
            cubemask=zeros(size(cube));
            cubemask(cube>maskat)=1;
            im.data=flip(permute(cubemask,[2 1 3]),3);
            if dilation>0, im.data=imdilate(im.data,se); end
            wfits(im,['masks/from_matlab/mask_' galname '_' vmethod '_v' int2str(version) '.fits']);
        case 'flat'
            cube=repmat(comapconv,[1 1 length(vel)]);
            cubemask=zeros(size(cube));
            cubemask((cube>maskat)&(abs(vv)<vfac*vmaxg))=1;
            im.data=flip(permute(cubemask,[2 1 3]),3);
            if dilation>0, im.data=imdilate(im.data,se); end
            wfits(im,['masks/from_matlab/mask_' galname '_' vmethod '_v' int2str(version) '.fits']);
        case 'block'
            cubemask=zeros(size(vv));
            cubemask((abs(vv)<vfac*vmaxg)&(rr25<1))=1;
            im.data=flip(permute(cubemask,[2 1 3]),3);
            if dilation>0, im.data=imdilate(im.data,se); end
            wfits(im,['masks/from_matlab/mask_' galname '_' vmethod '_v' int2str(version) '.fits']);
        case 'rotnoHa'
            ellipsemask=zeros(size(vv));
            ellipsemask((abs(vv)<vfac*vmaxg)&(rr25<1))=1;
            dv=(vv-vcen);
            cube=ellipsemask.*exp(-dv.^2./(2*(vwid/2.35).^2));
            cubemask=zeros(size(cube));
            cubemask(cube>maskat)=1;
            %im.data=permute(cubemask,[2 1 3]);
            im.data=flip(permute(cubemask,[2 1 3]),3);
            wfits(im,['masks/from_matlab/mask_' galname '_' vmethod '+_v' int2str(version) '.fits']);
            vcen=mkvelfield(r,(r.crval(1)-RA*15)*cosd(r.crval(2))*3600,(r.crval(2)-DEC)*3600,PA+180,INC,R,velR);
            dv=(vv-vcen);
            cube=ellipsemask.*exp(-dv.^2./(2*(vwid/2.35).^2));
            cubemask=zeros(size(cube));
            cubemask(cube>maskat)=1;
            im.data=flip(permute(cubemask,[2 1 3]),3);
            if dilation>0, im.data=imdilate(im.data,se); end
            wfits(im,['masks/from_matlab/mask_' galname '_' vmethod '-_v' int2str(version) '.fits']);
        otherwise
            disp('Unknown method. Quitting.');
            return
    end
    
    figure(3);
    clf
    imagesc(ra,de,comapconv); axis image; set(gca,'xdir','rev','ydir','nor');
    xlabel('RA (arcsec)');
    ylabel('DE (arcsec)');
    colormap parula
    cb=colorbar; cb.Label.String='T_{MB} dV (K km/s)';
    title('Expected CO map from Ha emission')
    
    figure(2);
    clf
    imagesc(ra,de,vcen); axis image; set(gca,'xdir','rev','ydir','nor');
    xlabel('RA (arcsec)');
    ylabel('DE (arcsec)');
    colormap jet
    caxis(vmaxg*[-2 2]);
    cb=colorbar; cb.Label.String='Vcen (km/s)';
    title('Velocity map of Ha')
    drawnow
    
    if ngal==1 % do this only when using it for one galaxy, otherwise it takes too long
        figure(1);
        clf
        for i=1:length(vel)
            imagesc(((1:im.numpt(1))-im.crpix(1))*im.cdelt(1)*3600,((1:im.numpt(2))-im.crpix(2))*im.cdelt(2)*3600,im.data(:,:,i)');
            axis image
            set(gca,'xdir','rev','ydir','nor');
            caxis([0 1]);
            xlabel('RA (arcsec)');
            ylabel('DE (arcsec)');
            cvel=((1:im.numpt(3))-im.crpix(3))*im.cdelt(3)+im.crval(3);
            ht=text(30,30,sprintf('%0.0f\n',cvel(i)),'color','w','fontsize',24,'fontwei','bold');
            drawnow;
        end
    end
end
