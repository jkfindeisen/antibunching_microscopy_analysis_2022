classdef progressbar < handle
    properties (Access=private)
        h
    end
    methods
        function obj = progressbar(name)
            if nargin == 0
                name = '';
            end
            obj.h = waitbar(0, name, 'Visible', 'off');
        end
        
        function set( obj, x)
            if ishandle( obj.h )
                waitbar(x, obj.h);
                obj.h.Visible = 'on';               
            end
        end
        
        function show( obj )
            if ishandle( obj. h)
                obj.h.Visible = 'on';
            end
        end        
        
        function hide( obj )
            if ishandle( obj. h)
                obj.h.Visible = 'off';
            end
        end
            
        function close( obj )
            if ishandle( obj.h )
                close( obj.h );
            end
        end
        
        function delete( obj )
            delete( obj.h );
        end
    end
end
    
