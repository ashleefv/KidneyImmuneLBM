test_i = [29,32,31,34,2];
load('data/listOfSensSortIdx.mat','listOfSensSortIdx');

% SpeciesIdx = 2:35;
% % non-sensitive parameters list
% NonSensIdx = setdiff(SpeciesIdx, SortIdx);
for i = 1:length(test_i)
    SortIdx = listOfSensSortIdx{1};
    g = test_i(i) == SortIdx;
    if sum(g)>0
        [num2str(test_i(i)) ' is sensitive for Number']
    end

    SortIdx = listOfSensSortIdx{3};
    g = test_i(i) == SortIdx;
    if sum(g)>0
        [num2str(test_i(i)) ' is sensitive for Diameter']
    end
end

    %% Number
    % SortIdx = listOfSensSortIdx{1};
    % SpeciesIdx = 2:35;
    % % non-sensitive parameters list
    % NonSensIdx = setdiff(SpeciesIdx, SortIdx);
    % g = test_i(i) == SortIdx;

    test_i = 31;%[31, 34, 2];
    i = 1;
    SortIdx = listOfSensSortIdx{1};
    g = test_i(i) == SortIdx;
    whichsubplot = 221;

        fig1 = openfig('fig7_EF.fig');
        subplot(whichsubplot)
        ax1 = gca; 
        h1 = findobj(ax1,'Type','line')
        spot = find(g>0)
            whichCurve=length(SortIdx)+1-spot
        
        addplot(h1,whichCurve,test_i(i),whichsubplot,ax1)
        
        fig2 = openfig('fig7.fig');
        subplot(whichsubplot)
        ax2 = gca; 
        h2 = findobj(ax2,'Type','line')        
        % No inhibition
        whichCurve = length(SortIdx)+1
        addplot(h2,whichCurve,test_i(i),whichsubplot,ax2)
        spot = find(g>0)
                if test_i == 34;
            %MLCP
            whichCurve=length(SortIdx)+2-spot
        else
            whichCurve=length(SortIdx)+1-spot
        end
        addplot(h2,whichCurve,test_i(i),whichsubplot,ax2)

        fig3 = openfig('figH.fig');
        subplot(whichsubplot)
        ax3 = gca; 
        h3 = findobj(ax3,'Type','line')        
        spot = find(g>0)
                if test_i == 34;
            %MLCP
            whichCurve=length(SortIdx)+2-spot
        else
            whichCurve=length(SortIdx)+1-spot
        end
        addplot(h3,whichCurve,test_i(i),whichsubplot,ax3)   

        fig4 = openfig('figJ.fig');
        subplot(whichsubplot)
        ax4 = gca; 
        h4 = findobj(ax4,'Type','line')        
        spot = find(g>0)
        if test_i == 34;
            %MLCP
            whichCurve=length(SortIdx)+2-spot
        else
            whichCurve=length(SortIdx)+1-spot
        end
        addplot(h4,whichCurve,test_i(i),whichsubplot,ax4)      
    % else
    %     [num2str(i) ' is non-sensitive for Number']
    %     fig1 = openfig('figG.fig');
    %     hold on
    %     subplot(whichsubplot)
    %     ax1 = gca; 
    %     h = findobj(ax1,'Type','line')
    %     % No inhibition
    %     whichCurve = length(NonSensIdx)+1
    %     addplot(h,whichCurve,test_i(i),whichsubplot,ax1)
    %     gnon = test_i(i) == NonSensIdx
    %     spot = find(gnon>0)
    %     whichCurve=length(NonSensIdx)+1-spot
    %     addplot(h,whichCurve,test_i(i),whichsubplot,ax1)
    %     close(fig1)
    %     clear h
    % end
    xlabel('Time (weeks)'); ylabel('Fenestration Number');legend('-DynamicLegend')
    %  xlim([0 20])


    test_i = [31];
    i = 1;
    SortIdx = listOfSensSortIdx{3};
    g = test_i(i) == SortIdx;
    whichsubplot = 223;

        fig1 = openfig('fig7_EF.fig');
        subplot(whichsubplot)
        ax1 = gca; 
        h1 = findobj(ax1,'Type','line')
        spot = find(g>0)
        whichCurve=length(SortIdx)+1-spot
        addplot(h1,whichCurve,test_i(i),whichsubplot,ax1)
        
        fig2 = openfig('fig7.fig');
        subplot(whichsubplot)
        ax2 = gca; 
        h2 = findobj(ax2,'Type','line')        
        % No inhibition
        whichCurve = length(SortIdx)+1
        addplot(h2,whichCurve,test_i(i),whichsubplot,ax2)
        spot = find(g>0)
        whichCurve=length(SortIdx)+1-spot
        addplot(h2,whichCurve,test_i(i),whichsubplot,ax2)

        fig3 = openfig('figH.fig');
        subplot(whichsubplot)
        ax3 = gca; 
        h3 = findobj(ax3,'Type','line')        
        spot = find(g>0)
        whichCurve=length(SortIdx)+1-spot
        addplot(h3,whichCurve,test_i(i),whichsubplot,ax3)   

        fig4 = openfig('figJ.fig');
        subplot(whichsubplot)
        ax4 = gca; 
        h4 = findobj(ax4,'Type','line')        
        spot = find(g>0)
        whichCurve=length(SortIdx)+1-spot
        addplot(h4,whichCurve,test_i(i),whichsubplot,ax4)   

   xlabel('Time (weeks)'); ylabel('Fenestration Diameter (nm)');ylim([45 80]); legend('-DynamicLegend')

   labelstring = {'A','', 'B'};
   for v = 1:2:3
      subplot(2,2,v)
         hold on
         text(-0.15, 1.1, labelstring(v)', 'Units', 'normalized', 'FontWeight', 'bold','FontSize',8)
         set(gca,'FontName','Arial','FontSize',6)
   end
   figname = 'Fig50';
widthInches = 6.5;
    heightInches = 4.23;
    run('ScriptForExportingImages.m')
end



