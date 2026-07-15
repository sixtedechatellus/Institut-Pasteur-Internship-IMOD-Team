for d in /data/ValidationEcoli/*/; do
    if [ -f "${d}decoys_generated" ]; then
        mv "${d}decoys_generated" "${d}decoys_final.ism"
        echo "✅ Renommé : ${d}decoys_generated → decoys_final.ism"
    fi
done
