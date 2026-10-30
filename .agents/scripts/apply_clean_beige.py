for path in ['web/src/components/ExperimentsView.jsx', 'web/src/components/HackathonView.jsx', 'web/src/components/StudioView.jsx']:
    with open(path, 'r', encoding='utf-8') as f:
        text = f.read()
    text = text.replace("border: '1px solid var(--border-subtle)',", "")
    text = text.replace("border: '1px solid var(--border-subtle)'", "")
    text = text.replace("border: '1px solid currentColor',", "")
    text = text.replace("border: '1px solid rgba(180, 83, 9, 0.25)',", "")
    text = text.replace("border: '1px solid rgba(52, 211, 153, 0.35)',", "")
    text = text.replace("background: 'var(--bg-elevated)'", "background: 'var(--bg-dark)'")
    with open(path, 'w', encoding='utf-8') as f:
        f.write(text)
print("Updated ExperimentsView, HackathonView, and StudioView to tonal surfaces.")
