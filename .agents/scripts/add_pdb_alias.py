with open('server.py', 'r', encoding='utf-8') as f:
    text = f.read()

target = '@app.get("/api/structure/{target_id}")\ndef get_structure_pdb'
replacement = '@app.get("/api/structure/{target_id}")\n@app.get("/api/campaign/pdb/{target_id}")\ndef get_structure_pdb'

if target in text:
    text = text.replace(target, replacement, 1)
    with open('server.py', 'w', encoding='utf-8') as f:
        f.write(text)
    print('Updated server.py successfully')
else:
    # try CRLF
    target_crlf = '@app.get("/api/structure/{target_id}")\r\ndef get_structure_pdb'
    replacement_crlf = '@app.get("/api/structure/{target_id}")\r\n@app.get("/api/campaign/pdb/{target_id}")\r\ndef get_structure_pdb'
    if target_crlf in text:
        text = text.replace(target_crlf, replacement_crlf, 1)
        with open('server.py', 'w', encoding='utf-8') as f:
            f.write(text)
        print('Updated server.py with CRLF successfully')
    else:
        print('Target not found')

