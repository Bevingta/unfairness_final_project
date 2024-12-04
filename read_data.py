import os
import json
import xml.dom.minidom

def get_downloaded_names():
    """Returns a list of file names from downloaded.txt"""
    try:
        with open('downloaded.txt', 'r') as f:
            return f.readlines()
    except FileNotFoundError:
        print("Error: downloaded.txt file not found")
        return []

def read_file_data(file_name):
    """
    Reads and displays the contents of files from GDC API
    """
    file_name = file_name.strip()
    full_path = os.path.join("Downloads", file_name)
    
    print(f"\nReading: {full_path}")
    
    if not os.path.exists(full_path):
        print(f"Error: File {full_path} does not exist")
        return
    
    if os.path.getsize(full_path) == 0:
        print(f"Error: {file_name} is empty (0 bytes)")
        return
    
    try:
        with open(full_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
            # Try to parse as JSON first
            try:
                json_data = json.loads(content)
                print("File contains JSON data:")
                print(json.dumps(json_data, indent=2)[:1000])  # Print first 1000 chars
                if len(content) > 1000:
                    print("... (output truncated)")
                return
            except json.JSONDecodeError:
                pass
            
            # Try to parse as XML
            try:
                xml_dom = xml.dom.minidom.parseString(content)
                pretty_xml = xml_dom.toprettyxml()
                print("File contains XML data:")
                print(pretty_xml[:1000])  # Print first 1000 chars
                if len(pretty_xml) > 1000:
                    print("... (output truncated)")
                return
            except xml.parsers.expat.ExpatError:
                pass
            
            # If not JSON or XML, treat as plain text
            print("File contains plain text:")
            print(content[:1000])  # Print first 1000 chars
            if len(content) > 1000:
                print("... (output truncated)")
            
    except Exception as e:
        print(f"Error reading file: {str(e)}")

def main():
    downloaded_files = get_downloaded_names()
    if not downloaded_files:
        print("No files found to process")
        return
        
    for file_name in downloaded_files:
        read_file_data(file_name)
        print("\n" + "-" * 80)

if __name__ == "__main__":
    main()
