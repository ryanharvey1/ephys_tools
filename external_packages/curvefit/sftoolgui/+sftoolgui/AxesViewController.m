classdef AxesViewController < handle
    %AxesViewController controls the axes view
    %
    %   AxesViewController(AXES, AXESVIEWMODEL, SHOW2DVIEW) controls the
    %   AXES View property. It disables Rotate3D mode if the view is 2D and
    %   enables it otherwise.
    
    %   Copyright 2011 The MathWorks, Inc.
    %   $Revision: 1.1.6.2.2.1 $    $Date: 2011/07/18 00:30:46 $
    
    properties( SetAccess = 'public', GetAccess = 'public', Dependent)
        % View2D is true if the panel is currently displaying curve
        % data and false otherwise.
        View2D = false;
    end
    
    properties(SetAccess = 'private', GetAccess = 'private')
        % HAxes is the handle to the axes being controlled
        HAxes ;
        % PrivateView2D is true if the panel is currently displaying curve
        % data and false otherwise. It is associated with the public View2D
        % dependent property.
        PrivateView2D = false;
        % AxesViewModel is the sftoolgui.AxesViewModel that maintains the
        % view angle for the axes under control.
        AxesViewModel;
        % ModelViewChangedListener listens for changes on the axes View
        % property not initiated in this class.
        ModelViewChangedListener;
        % ViewPropertyListener is a listener on the axes View property
        ViewPropertyListener;
    end
    
    methods
        function this = AxesViewController(axes, axesViewModel, show2DView)
            % AxesViewController is the constructor for the Curve Fitting
            % Tool Axes View Property Controller
            %
            % this = AxesViewController(AXES, AXESVIEWMODEL, SHOW2DVIEW)
            % creates an AxesViewController that controls AXES.
            % AXESVIEWMODEL is the an sftoolgui.AxesViewModel. SHOW2DVIEW
            % is true if the axes should show 2D data and false otherwise.
            
            % Make sure first input is an axes
            if ~ishghandle(axes,'axes')
                error(message('curvefit:AxesViewController:InvalidInput'));
            end
            
            % Make sure second argument is an AxesViewModel
            if ~isa(axesViewModel, 'sftoolgui.AxesViewModel')
                error(message('curvefit:AxesViewController:InvalidInputAxesViewModel'));
            end
            
            % Make sure third argument is a logical
            if ~islogical(show2DView)
                error(message('curvefit:AxesViewController:InvalidInputShow2DView'));
            end
            
            this.AxesViewModel = axesViewModel;
            this.HAxes = axes;
            
            % Listen to the axes View property changing from other actions
            % outside this class, for instance when users rotate the axes
            this.ViewPropertyListener = curvefit.proplistener( this.HAxes, 'View', @(s, e) this.axesViewChanged );
            
            % Listen to AxesViewModel 'ThreeDViewAngleChanged' event
            this.ModelViewChangedListener = event.listener(this.AxesViewModel, 'ThreeDViewAngleChanged', @(s, e) this.axesViewModelThreeDViewAngleChanged(  e ) );
            
            % Initialize the axes and default view
            initializeViews(this, show2DView);
        end
        
        function view2D = get.View2D(this)
            % Return the value of the associated private property.
            view2D = this.PrivateView2D;
        end
        
        function set.View2D(this, show2DView)
            %In addition to setting the associated PrivateView2D property,
            %this set method handles switching between 2D and 3D views. It
            %sets the axes View property, updates plotview information and
            %enables or disables rotation.
            if show2DView
                defaultViewAngle = sftoolgui.util.DefaultViewAngle.TwoD;
                if ~this.PrivateView2D
                    setViewQuietly(this, defaultViewAngle);
                end
            else
                defaultViewAngle = sftoolgui.util.DefaultViewAngle.ThreeD;
                if this.PrivateView2D
                    % Switching from curve to non curve, restore the
                    % previous surface view.
                    setViewQuietly(this, this.AxesViewModel.ThreeDViewAngle);
                end
            end
            % Update the plotview information
            setDefaultView(this, defaultViewAngle);
            % Update PrivateView2D to reflect state
            this.PrivateView2D = show2DView;
            % Enable/disable Rotation based on view dimension.
            behavior = hggetbehavior( this.HAxes, 'Rotate3d' );
            set( behavior, 'Enable', ~this.PrivateView2D );
        end
    end
    
    methods(Access = private)
        function initializeViews(this, show2DView)
            % initializeViews initialize the views.
            %
            %    initializeView(THIS, SHOW2DVIEW) will set THIS.HAXES view
            %    property to the curve values if SHOW2DVIEW is true or to
            %    surface values otherwise.
            if show2DView
                viewAngle = sftoolgui.util.DefaultViewAngle.TwoD;
                setViews(this, viewAngle, viewAngle);
            else
                setViews(this,...
                    this.AxesViewModel.ThreeDViewAngle, ...
                    sftoolgui.util.DefaultViewAngle.ThreeD);
            end
            this.View2D = show2DView;
        end
        
        function setViews(this, viewAngle, defaultViewAngle)
            % setViews will set the axes' View Property without firing an event
            % and will set the Rotate3D's default view.
            setViewQuietly(this, viewAngle);
            setDefaultView(this, defaultViewAngle);
        end
        
        function setViewQuietly( this, viewAngle )
            % setViewQuietly quietly sets the axes View property.
            %
            %    setViewQuietly(THIS, VIEWANGLE) disables the listener on
            %    THIS.HAXES View property, sets the View property and then
            %    (re)enables the listener.
            curvefit.setListenerEnabled(this.ViewPropertyListener, false);
            set(this.HAxes, 'View', viewAngle);
            curvefit.setListenerEnabled(this.ViewPropertyListener, true);
        end
        
        function setDefaultView( this, defaultView )
            %setDefaultView sets Rotate3D's default view.
            %
            %   setDefaultView(THIS, DEFAULTVIEW) sets THIS.HAXES plotview
            %   information's View property to DEFAULTVIEW.
            %
            %   cftool wants to control rotate3D mode's "Reset to Original
            %   View" functionality. Specifically, when users choose that
            %   option, the View property should be changed to the plot
            %   panel's default view, which may change depending on
            %   selected data. Rotate3D's view information, which can be
            %   set with a call to 'resetplotview(axes,
            %   'SaveCurrentView')', needs to be updated with the default
            %   view to achieve this result.
            
            % Get the user's View.
            currentView = get(this.HAxes, 'View');
            
            setViewQuietly(this, defaultView);
            
            % Save the information.
            resetplotview(this.HAxes, 'SaveCurrentView' );
            
            setViewQuietly(this, currentView);
        end
        function axesViewModelThreeDViewAngleChanged(this, evt)
            % axesViewModelThreeDViewAngleChanged is the callback for sftoolgui.AxesViewModel ThreeDViewAngleChanged event.
            %
            % It updates the View property. If View2D is true, it sets the axes
            % to the default 2D angle. Otherwise, it sets the View property the
            % new 3D angle.
            if this.View2D
                viewAngle = this.AxesViewModel.Default2DViewAngle;
            else
                viewAngle = evt.AxesViewModel.ThreeDViewAngle;
            end
            setViewQuietly(this, viewAngle);
        end
        
        function axesViewChanged(this)
            % axesViewChanged is the callback when the View property has been
            % changed from actions outside this class (which only happens with
            % 3D views).
            %
            % It sets the AxesViewModels "ThreeDViewAngle" property.
            this.AxesViewModel.ThreeDViewAngle = get(this.HAxes, 'View');
        end
    end
end

