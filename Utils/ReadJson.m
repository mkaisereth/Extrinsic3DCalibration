classdef ReadJson
    % Class to Read JSON
    
    methods(Static)
        
        function fileContent = ReadTest()
            % Read a test json
            % return the file content as decoded json
            
            fname = 'Test.json';
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            % decode json
            fileContent = jsondecode(str);
        end
        
        function fileContent = Read(filePath)
            % Read a json file from filePath
            % filePath - the path to the json file
            % returns the file content as decoded json
            
            fid = fopen(char(filePath));
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            
            str = strrep(str,"\\","\");
            str = strrep(str,"\","\\");
            
            % decode json
            fileContent = jsondecode(str);
        end
        
        function Write(filePath, jsonObj)
            % Write a json file from an object to filePath
            % filePath - the path to the json file
            % jsonObj - the json object
            
            fileContent = jsonencode(jsonObj);
            % add a return character after all commas:
            fileContent = strrep(fileContent, ',', ',\n');
            % add a return character after curly brackets:
            fileContent = strrep(fileContent, '{', '{\n');
            fileContent = strrep(fileContent, '}', '\n}');
            
            fid = fopen(char(filePath), 'w');
            fprintf(fid, fileContent);
            fclose(fid);
        end
    end
end