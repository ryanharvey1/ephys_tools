classdef LimitSpinnersDialog < handle
    %Surface Fitting Tool Limit Spinners Dialog
    
    %   Copyright 2008-2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.6.2.1 $    $Date: 2011/07/18 00:30:53 $
    
    properties (SetAccess = 'private', GetAccess = 'private')
        HFitFigure ;
        JavaDialog ;
        InternalListeners ;
        FitdevListeners ;
        AxesLimitListeners ;
    end
    
    properties(Dependent = true, SetAccess = 'private')
        XLim ;
        SurfaceYLim ;
        CurveYLim ;
        ResidYLim ;
        MainZLim ;
        ResidZLim ;
    end
    
    methods
        function this = LimitSpinnersDialog(hFitFigure)
            this.HFitFigure = hFitFigure;
            
            this.JavaDialog = javaObjectEDT('com.mathworks.toolbox.curvefit.surfacefitting.AxesLimitsDialog');
            
            this.InternalListeners = [
                addlistener( this.JavaDialog, 'resetLimits', @(s, e) this.resetDefaultAxisLimits() );
                addlistener( this.JavaDialog, 'limitsChanged', @(s, e) this.axesLimitsDialogAction( e ));
                ];
            
            this.FitdevListeners = [
                event.listener(this.HFitFigure.HFitdev, 'FittingDataUpdated', @(s, e) this.fittingDataUpdated(  ) );
                event.listener(this.HFitFigure.HFitdev, 'FitNameUpdated', @(s, e) this.fitNameUpdated( ) );
                event.listener(this.HFitFigure.HFitdev, 'DimensionChanged', @(s, e) this.dimensionChangedAction() );
                ];
            
            this.AxesLimitListeners = event.listener(this.HFitFigure.AxesViewModel, ...
                'LimitsChanged', @(s, e) this.updateAllLimitSpinners( e ) );
            
            % Set the title
            fitNameUpdated(this);
            
            % Set variable names
            fittingDataUpdated(this);
            
            % Make sure the correct view is showing
            dimensionChangedAction(this);
            
            % Set spinner values.
            updateAllLimitSpinners(this, sftoolgui.AxesViewEventData(this.HFitFigure.AxesViewModel));
        end
        
        function show(this)
            % show makes the dialog visible
            javaMethodEDT('showDialog', this.JavaDialog);
        end
        
        function updateAllLimitSpinners(this, evt)
            % updateAllLimitSpinners updates all the java spinner widgets
            if isCurveDataSpecified(this.HFitFigure.HFitdev.FittingData)
                setValueAndStep(this, 'CurveYLim', evt.AxesViewModel.ResponseLimits);
            end
            
            setValueAndStep(this, 'XLim', evt.AxesViewModel.XInputLimits);
            setValueAndStep(this, 'SurfaceYLim', evt.AxesViewModel.YInputLimits);
            setValueAndStep(this, 'ResidYLim', evt.AxesViewModel.ResidualLimits);
            setValueAndStep(this, 'MainZLim', evt.AxesViewModel.ResponseLimits);
            setValueAndStep(this, 'ResidZLim', evt.AxesViewModel.ResidualLimits);
        end
        
        function delete(this)
            % This is the class's delete method
            javaMethodEDT('dispose', this.JavaDialog);
            deleteListeners(this);
        end
        
        function fittingDataUpdated(this)
            % fittingDataUpdated sets the spinner labels
            [xname, yname, zname] = getNames(this.HFitFigure.HFitdev.FittingData);
            javaMethodEDT('setLabels', this.JavaDialog, xname, yname, zname)
        end
        
        function dimensionChangedAction(this)
            % dimensionChangedAction lets the java dialog know that the
            % dimension has changed.
            if isCurveDataSpecified(this.HFitFigure.HFitdev.FittingData)
                javaMethodEDT('showDimensionView', this.JavaDialog, 2);
            else
                javaMethodEDT('showDimensionView', this.JavaDialog, 3);
            end
        end
        
        function fitNameUpdated(this)
            % fitNameUpdated updates the fit name on the java dialog title
            title = xlate('Axes Limits - ');
            title = sprintf('%s%s', title, this.HFitFigure.HFitdev.FitName);
            javaMethodEDT('setTitle', this.JavaDialog, title);
        end
        
        function lim = get.XLim(this)
            lim = this.HFitFigure.AxesViewModel.XInputLimits;
        end
        
        function set.XLim(this, lim)
            this.HFitFigure.AxesViewModel.setLimit('XInputLimits', lim);
        end
        
        function lim = get.SurfaceYLim(this)
            lim = this.HFitFigure.AxesViewModel.YInputLimits;
        end
        
        function set.SurfaceYLim(this, lim)
            this.HFitFigure.AxesViewModel.setLimit('YInputLimits', lim);
        end
        
        function lim = get.CurveYLim(this)
            lim = this.HFitFigure.AxesViewModel.ResponseLimits;
        end
        
        function set.CurveYLim(this, lim)
            this.HFitFigure.AxesViewModel.setLimit('ResponseLimits', lim);
        end
        
        function lim = get.ResidYLim(this)
            lim = this.HFitFigure.AxesViewModel.ResidualLimits;
        end
        
        function set.ResidYLim(this, lim)
            this.HFitFigure.AxesViewModel.setLimit('ResidualLimits', lim);
        end
        
        function lim = get.MainZLim(this)
            lim = this.HFitFigure.AxesViewModel.ResponseLimits;
        end
        
        function set.MainZLim(this, lim)
            this.HFitFigure.AxesViewModel.setLimit('ResponseLimits', lim);
        end
        
        function lim = get.ResidZLim(this)
            lim = this.HFitFigure.AxesViewModel.ResidualLimits;
        end
        
        function set.ResidZLim(this, lim)
            this.HFitFigure.AxesViewModel.setLimit('ResidualLimits', lim);
        end
        
        function deleteListeners(this)
            % deleteListeners deletes the java dialog and Fitdev listeners
            delete(this.FitdevListeners);
            delete(this.InternalListeners);
        end
    end
    
    methods(Access = 'private')
        
        function resetDefaultAxisLimits(this)
            % resetDefaultAxisLimits calls plotData which will set limits
            % based on the data.
            plotData(this.HFitFigure, true);
        end
        
        function axesLimitsDialogAction(this, e)
            % axesLimitsDialogAction sets the corresponding property of the
            % spinner changed.
            property = char(e.getLimitProperty());
            this.(property) = [e.getMinValue() e.getMaxValue()];
        end
        
        function setValueAndStep(this, spinner, limit)
            % setValueAndStep sets the java spinners value and step
            javaMethodEDT('setValueAndStep', this.JavaDialog, spinner, ...
                limit(1), limit(2), iCalculateStep(limit));
        end
    end
end

function step = iCalculateStep(lim)
% iCalculateStep calculates the step.
dx = (lim(2)-lim(1))/100;
pwr = floor(log10(dx));
step = 10^pwr * (round(dx/(10^pwr)));
end
