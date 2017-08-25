function [abaPath,cmissPath,elemSize,dispOfInterest] = ...
    case_selector(LOADING,REFINEMENT,DIMENSION,INTERPOLATION,CONTROL)


abaPrefix = 'reference/abaqus/';
ironPrefix = 'reference/iron/';

switch DIMENSION
    case '3D'
        switch LOADING
            case 'SHEAR'
                dispOfInterest = 1; %x
                switch REFINEMENT
                    case 'COARSE'
                        if strcmp(CONTROL, 'DISPLACEMENT')
                            abaPath = [abaPrefix '3D_SHEAR_ELASTIC_elem_20_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n08x06x06_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_disp.txt'];                            
                        elseif strcmp(CONTROL, 'FORCE')
                            abaPath = [abaPrefix '3D_SHEAR_FORCE_ELASTIC_elem_20_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n08x06x06_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_force.txt'];   
                        end
                        elemSize = 20;
                    case 'MEDIUM'
                            if strcmp(CONTROL, 'DISPLACEMENT')                        
                            abaPath = [abaPrefix '3D_SHEAR_ELASTIC_elem_10_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n16x12x12_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_disp.txt'];                            
                        elseif strcmp(CONTROL, 'FORCE')                        
                            abaPath = [abaPrefix '3D_SHEAR_FORCE_ELASTIC_elem_10_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n16x12x12_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_force.txt'];                            
                        end                        
                        elemSize = 10;
                    case 'FINE'
                        if strcmp(CONTROL, 'DISPLACEMENT')                        
                            abaPath = [abaPrefix '3D_SHEAR_ELASTIC_elem_5_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n32x24x24_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_disp.txt'];                            
                        elseif strcmp(CONTROL, 'FORCE')                        
                            abaPath = [abaPrefix '3D_SHEAR_FORCE_ELASTIC_elem_5_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n32x24x24_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_force.txt'];
                        end                        
                        elemSize = 5;
                end
            case 'UNIAX'
                dispOfInterest = 2; %x
                switch REFINEMENT
                    case 'COARSE'
                        if strcmp(CONTROL, 'DISPLACEMENT')                        
                            abaPath = [abaPrefix '3D_UNIAX_ELASTIC_elem_20_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n08x06x06_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_disp.txt'];                            
                        elseif strcmp(CONTROL, 'FORCE')                        
                            abaPath = [abaPrefix '3D_UNIAX_FORCE_ELASTIC_elem_20_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n08x06x06_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_force.txt'];
                        end                        
                        elemSize = 20;
                    case 'MEDIUM'
                        if strcmp(CONTROL, 'DISPLACEMENT')                        
                            abaPath = [abaPrefix '3D_UNIAX_ELASTIC_elem_10_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n16x12x12_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_disp.txt'];                            
                        elseif strcmp(CONTROL, 'FORCE')                        
                            abaPath = [abaPrefix '3D_UNIAX_FORCE_ELASTIC_elem_10_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n16x12x12_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_force.txt'];
                        end                        
                        elemSize = 10;
                    case 'FINE'
                        if strcmp(CONTROL, 'DISPLACEMENT')                        
                            abaPath = [abaPrefix '3D_UNIAX_ELASTIC_elem_5_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n32x24x24_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_disp.txt'];                            
                        elseif strcmp(CONTROL, 'FORCE')                        
                            abaPath = [abaPrefix '3D_UNIAX_FORCE_ELASTIC_elem_5_160x120x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x120_n32x24x24_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_force.txt'];                            
                        end                        
                        elemSize = 5;                        
                end
        end
    case '2D'
        switch LOADING
            case 'SHEAR'
                dispOfInterest = 1; %x
                switch REFINEMENT
                    case 'COARSE'
                        if strcmp(CONTROL, 'DISPLACEMENT')                        
                            abaPath = [abaPrefix '2D_SHEAR_ELASTIC_elemSize_20_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x000_n08x06x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_disp.txt'];                            
                        elseif strcmp(CONTROL, 'FORCE')                        
                            abaPath = [abaPrefix '2D_SHEAR_FORCE_ELASTIC_elemSize_20_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x000_n08x06x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_force.txt'];     
                            
                        end                        
                        elemSize = 20;
                    case 'MEDIUM'
                        if strcmp(CONTROL, 'DISPLACEMENT')                        
                            abaPath = [abaPrefix '2D_SHEAR_ELASTIC_elemSize_10_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x000_n16x12x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_disp.txt'];                            
                        elseif strcmp(CONTROL, 'FORCE')                        
                            abaPath = [abaPrefix '2D_SHEAR_FORCE_ELASTIC_elemSize_10_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x000_n16x12x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_force.txt'];   
                            
                        end                        
                        elemSize = 10;
                    case 'FINE'
                        if strcmp(CONTROL, 'DISPLACEMENT')                        
                            abaPath = [abaPrefix '2D_SHEAR_ELASTIC_elemSize_5_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x000_n32x24x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_disp.txt'];                            
                        elseif strcmp(CONTROL, 'FORCE')                        
                            abaPath = [abaPrefix '2D_SHEAR_FORCE_ELASTIC_elemSize_5_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x000_n32x24x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_shear_force.txt'];
                            
                        end                        
                        elemSize = 5;
                end
            case 'UNIAX'
                dispOfInterest = 2; %x
                switch REFINEMENT
                    case 'COARSE'
                        if strcmp(CONTROL, 'DISPLACEMENT')                        
                            abaPath = [abaPrefix '2D_UNIAX_ELASTIC_elem_20_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x000_n08x06x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_disp.txt'];                            
                        elseif strcmp(CONTROL, 'FORCE') 
                            abaPath = [abaPrefix '2D_UNIAX_FORCE_ELASTIC_elem_20_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = [ 'l160x120x000_n08x06x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_force.txt'];                               
                        end                        
                        elemSize = 20;
                    case 'MEDIUM'
                        if strcmp(CONTROL, 'DISPLACEMENT')                        
                            abaPath = [abaPrefix '2D_UNIAX_ELASTIC_elem_10_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = ['l160x120x000_n16x12x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_disp.txt'];                            
                        elseif strcmp(CONTROL, 'FORCE')                        
                            abaPath = [abaPrefix '2D_UNIAX_FORCE_ELASTIC_elem_10_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = ['l160x120x000_n16x12x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_force.txt'];
                            
                        end                                                
                        elemSize = 10;
                    case 'FINE'
                        if strcmp(CONTROL, 'DISPLACEMENT')                        
                            abaPath = [abaPrefix '2D_UNIAX_ELASTIC_elem_5_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = ['l160x120x000_n32x24x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_disp.txt'];
                        elseif strcmp(CONTROL, 'FORCE')                        
                            abaPath = [abaPrefix '2D_UNIAX_FORCE_ELASTIC_elem_5_160x120mm_intp_' num2str(INTERPOLATION) '_DIRECT.txt'];
                            tmp = ['l160x120x000_n32x24x00_i' num2str(INTERPOLATION) '_s0'];
                            cmissPath = [ironPrefix tmp '/' tmp '_uniax_force.txt'];
                            
                        end                        
                        elemSize = 5;                        
                end               
        end        
end