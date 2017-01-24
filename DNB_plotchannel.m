function DNB_plotchannel(b_filter,a_filter,filt_type,ChannelNr)
%Plotting the signal before and after filtering for every channel
%b_filter = signal before filtering
%a_filter = signal after filtering
%filt_type = type of filtering. Options: LPF, HPF, all, MARA, noMARA
%           LPF = low pass filtering: LPF, gaussian or hrf filtering
%           HPF = high pass filtering: detrending, wavelet or DCT filtering
%           all = both (MARA), low and high pass filtering
%           MARA = Movement Artifact Removal Algorithm
%           noMARA = the DNB_MARA script was used over the signal, but it did not
%           change the signal (moving standard deviation signal does not reach
%           threshold).
%ChannelNr = Number of channel plotted. Can only be used when plotting separate channels; when all channels are entered in the function simultaneously, ChannelNr should be left open.
%
% Lars Didden - Donders Centre for Cognitive Neuroimaging
% Joost Wegman - Donders Centre for Cognitive Neuroimaging

global INFO;

isess=INFO.counter.iSess;

try 
    aa_filter=a_filter{isess};
catch
    aa_filter=a_filter;
end

if numel(b_filter) ~= numel(aa_filter)
    disp('The size of the signal before and after filtering is not equal.')
end

for iChan=1:size(b_filter,2);
    figure;
    if strcmp(filt_type,'noMARA')==1
        h1=plot(b_filter(:,iChan),'m');
        clear h1
    elseif strcmp(filt_type,'MARA')==1
        h1=plot(b_filter(:,iChan),'b'); hold on; h2=plot(aa_filter(:,iChan),'r');
        clear h1 h2
    elseif strcmp(filt_type,'LPF')==1
        diffba=b_filter(:,iChan)-aa_filter(:,iChan);
        h1=plot(diffba,'k'); hold on; h2=plot(b_filter(:,iChan),'b'); h3=plot(aa_filter(:,iChan),'r');
        set(h1,'Linewidth',0.5); set(h2,'Linewidth',2); set(h3,'Linewidth',2);
        clear h1 h2 h3 diffba
    elseif strcmp(filt_type,'HPF')==1
        h1=plot(b_filter(:,iChan),'b'); hold on; h2=plot(aa_filter(:,iChan),'r');
        set(h1,'Linewidth',2); set(h2,'Linewidth',2);
        clear h1 h2
    else
        h1=plot(b_filter(:,iChan),'b'); hold on; h2=plot(aa_filter(:,iChan),'r');
        set(h2,'Linewidth',2);
        clear h1 h2
    end
    
    if strcmp(INFO.model.hb,'HbO')==1; hb='Oxyhemoglobin';
    elseif strcmp(INFO.model.hb,'HbR')==1; hb='Deoxyhemoglobin';
    end
    
    if strcmp(INFO.SCI.check,'yes')==1 
        if nargin==4
            title(['\fontsize{15}',hb,' channel ',num2str(INFO.SCI.sessremchannel{end}(ChannelNr)),]);
            plotname=['channel ',num2str(INFO.SCI.sessremchannel{end}(ChannelNr)),];
        else
            title(['\fontsize{15}',hb,' channel ',num2str(INFO.SCI.sessremchannel{end}(iChan)),]);
            plotname=['channel ',num2str(INFO.SCI.sessremchannel{end}(iChan)),];
        end
    else
        if nargin==4
            title(['\fontsize{15}',hb,' channel ',num2str(ChannelNr),]);
            plotname=['channel ',num2str(ChannelNr),];
        else
            title(['\fontsize{15}',hb,' channel ',num2str(iChan),]);
            plotname=['channel ',num2str(iChan),];
        end
    end
    
    if strcmp(filt_type,'HPF')==1; legend(['Before ',INFO.filt.HPF,' filtering'],['After ',INFO.filt.HPF,' filtering']);
        plotroot=fullfile(INFO.file.figures{isess},INFO.file.name,'PP_afterHPF');
    elseif strcmp(filt_type,'LPF')==1; legend(['Difference between before and after ',INFO.filt.LPF,],['Before ',INFO.filt.LPF,' filtering'],['After ',INFO.filt.LPF,' filtering']);
        plotroot=fullfile(INFO.file.figures{isess},INFO.file.name,'PP_afterLPF');
    elseif strcmp(filt_type,'all')==1; legend('Before filtering','After filtering');
        plotroot=fullfile(INFO.file.figures{isess},INFO.file.name,'PP_before_after_filtering');
    elseif strcmp(filt_type,'MARA')==1; legend('Before MARA','After MARA');
        plotroot=fullfile(INFO.file.figures{isess},INFO.file.name,'MARA_figures');
    elseif strcmp(filt_type,'noMARA')==1; legend('Before and after MARA');
        plotroot=fullfile(INFO.file.figures{isess},INFO.file.name,'MARA_figures');
    end
    
    if ~exist(plotroot); mkdir(plotroot); end
    
    plotdir=fullfile(plotroot,plotname);
    saveas(gcf,plotdir,INFO.extension);
    close(gcf);
end



