%demander nom du fichier 1 � impoter [nameC1,filepathC1]=uigetfile('*.*','Select file cell 1');filename=fullfile(filepathC1,nameC1);data=importdata(nameC1,'\t',1);%Importer un fichier provenant de read (labview) %Cell 1%s�parer les diff�rentes variables pour cell 1.tps=data.data(:,1);redx=data.data(:,2);redy=data.data(:,3);greenx=data.data(:,4);greeny=data.data(:,5);spike=data.data(:,6); % col (:,6) donne les spikes du spike counter 1, col (:,7) ceux du spike counter 2angle=data.data(:,10);%avec Col 1=temps, col 2=x led rouge, col 3=y led rouge, col 4=x led verte,%col 5=y led verte, col 6=spike or no spyke, col 10 = angle en radians, col 11 =angle en degr�s.%Cell 1%elimine les points o� la diode rouge n'est pas d�tect�e pour cell 1bads=find(redx==0);greenx2=greenx;greenx2(bads)=[];greeny2=greeny;greeny2(bads)=[];spike2=spike;spike2(bads)=[];angle2=angle;angle2(bads)=[];tps2=tps;tps2(bads)=[];%Cell 1%elimine les points o� la diode verte n'est pas d�tect�ebads2=find(greenx2==0);greenx3=greenx2;greenx3(bads2)=[];greeny3=greeny2;greeny3(bads2)=[];spike3=spike2;spike3(bads2)=[];angle3=angle2;angle3(bads2)=[];tps3=tps2;tps3(bads2)=[];%Cell 1%partage le cercle en 60 parties (6� chacune) pour cell 1nBins=60;for i=1:nBins;% Nombre d'instance o� la t�te du rat s'est trouv� dans chacune de ces 45% possibles orientations.    ListOrientation=find(angle3>=((((2*pi)/nBins)/2)+(i-1)*2*pi/nBins) & angle3<((((2*pi/nBins)/2)+(i)*2*pi/nBins)));if length(ListOrientation)<1;%dans les cas ce nombre vaut 0, on attribue la valeur afin d'�viter une%division par 0    NombreOrientation(i)=1;% et on attribue automatiquement 0 au nombre de spikes pour cette% orientation.    NombreSpikesOrientation(i)=0;      elseNombreSpikesOrientation(i)=sum(spike3(ListOrientation));NombreOrientation(i)=length(ListOrientation); endend% transforme les valeurs (60i�me de sec) en nombre de spikes par seconde%Cell 1NombreSpikesOrientation2=NombreSpikesOrientation*60; % *60ms pour avoir le temps en s% Calcule la frequence de d�charge (spike/sec) pour chaque bin de 8�.BinsNbSpikes=NombreSpikesOrientation2./NombreOrientation;%trace du graphe du taux de d�charge/direction pour les cellules 1BinsAngle=[0+((2*pi/(nBins)/2)):2*pi/(nBins):(2*pi)-((2*pi/(nBins)/2))];% subplot(2,2,1),h=polarKB1Cell(BinsAngle,BinsNbSpikes,'k',1);% h=title(nameC1);% hold on;% Second trac� des point 45 et 1 afin de relier les 2 extr�mit�s de la% courbe pour cell 1 & 2BinsAngle2=[(2*pi)-(((2*pi)/nBins)/2):(((2*pi)/nBins))/2:(2*pi)+ ((2*pi)/nBins)/2];BinsNbSpikes2=[BinsNbSpikes(45) (BinsNbSpikes(45)+BinsNbSpikes(1))/2 BinsNbSpikes(1)];% subplot(2,2,1),hh=polarKB1Cell(BinsAngle2,BinsNbSpikes2,'k',1);%Affichage des taux de d�charge Cell 1MaxFiringRateC1=max(BinsNbSpikes);FRCell1 = sprintf('FR=%7.2f',MaxFiringRateC1);text(-50,60,FRCell1,'fontsize',6);hold on;BinsAngle3=BinsAngle*nBins;subplot(2,2,3:4),h2=plot(BinsAngle3,BinsNbSpikes);axis tight;xlim ([0 360]);hold on;%Trajet des animaux avec spikeshold on;subplot(2,2,2),e=plot(greenx3,greeny3);set(gca,'xtick',[],'ytick',[]);set(e,'Color','black','LineWidth',1);axis([0 255 0 255]);axis square;e1=title('trajet avec spikes');hold on;spike4=find(spike3>0);greenx4=greenx3(spike4);greeny4=greeny3(spike4);subplot(2,2,2),ee=plot(greenx4,greeny4,'.r');axis([0 255 0 255] );axis square;