classdef CTimeleft < handle
% How to use:
% t0 = CTimeleft(looplength,interval);
% for i = 1:interval:looplength
%   t0.timeleft();
%   ... do stuff ...
% end
    properties 
        t0
        charsToDelete = [];
        done
        total
        interval = 1;
    end
    
    methods
        function t = CTimeleft(total, interval)
            
            t.done = 0;
            t.total = total;
            if nargin>1
                t.interval = interval;
            else
                t.interval = ceil(total*t.interval/100);
            end
        end
        
        function [remaining, status_string] = timeleft(t)
            if t.done == 0
                t.t0 = tic;
            end
                              
            t.done = t.done + 1;
            
            elapsed = toc(t.t0);
            
            if t.done == 1 || mod(t.done,t.interval)==0 || t.done == t.total || nargout > 0

                % compute statistics
                avgtime = elapsed./t.done;
                remaining = (t.total-t.done)*avgtime;
                
                if avgtime < 1
                    ratestr = sprintf('- %.2f iter/s', 1./avgtime);
                else
                    ratestr = sprintf('- %.2f s/iter', avgtime);
                end
                
                if t.done == 1
                    remaining = -1;
                    ratestr = [];
                end
                
                timesofarstr  = t.sec2timestr(elapsed);
                timeleftstr = t.sec2timestr(remaining);
                
                %my beloved progress bar
                pbarw = 30;
                pbar = repmat('.',1,pbarw);
                pbarind = 1 + round(t.done/t.total*(pbarw));
                pbar(1:pbarind-1) = '=';
                if pbarind < pbarw
                    pbar(pbarind) = '>';
                else
                        0;
                    
                end
                pbar = ['[',pbar,']'];
                
                
                status_string = sprintf('%s %03d/%03d - %03d%%%% - %s|%s %s ',pbar,t.done,t.total,...
                    floor(t.done/t.total*100),timesofarstr,timeleftstr, ratestr);
                
                delstr = [];
                if ~isempty(t.charsToDelete)
                    delstr = repmat('\b',1,t.charsToDelete-1);
                end
           
                if nargout == 0
                    fprintf([delstr status_string]);
                end
                
                t.charsToDelete = numel(status_string);
            end
            
            if t.done == t.total && nargout == 0 
                fprintf('\n');
            end
        end
        
        function timestr = sec2timestr(~,sec)
            % Convert a time measurement from seconds into a human readable string.

            if sec < 0
                timestr = '???:???';
                return
            end

            % Convert seconds to other units
            d = floor(sec/86400); % Days
            sec = sec - d*86400;
            h = floor(sec/3600); % Hours
            sec = sec - h*3600;
            m = floor(sec/60); % Minutes
            sec = sec - m*60;
            s = floor(sec); % Seconds

            % Create time string
            if d > 0
                timestr = sprintf('%02dd:%02dh',d,h);
            elseif h > 0
                timestr = sprintf('%02dh:%02dm',h,m);
            else
                timestr = sprintf('%02dm:%02ds',m,s);
            end
        end
    end
end
   


