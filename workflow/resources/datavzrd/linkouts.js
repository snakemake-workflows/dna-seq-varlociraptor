function hgvsc_dropdown(value, row) {
  function varsome_link(row) {
    let ref = row["reference allele"];
    let alt = row["alternative allele"];
    if ( alt && ref) {
      let chr = row.chromosome;
      let pos = row.position;
      let build = ( row.build === "GRCh38" ) ? "hg38/" : "hg19/";
      const url = "https://varsome.com/variant/";
      let aa_regex = /^[ACGTacgt]+$/;
      if ( ref.match(aa_regex) && alt.match(aa_regex)) {
        let descriptor = `${chr}:${pos}:${ref}:${alt}`
        return `${url}${build}${descriptor}`;
      }
    }
    return ``
  }

  function genebe_link(row) {
    let ref = row["reference allele"];
    let alt = row["alternative allele"];
    if ( alt && ref) {
      let chr = row.chromosome;
      let pos = row.position;
      let build = ( row.build === "GRCh38" ) ? "hg38/" : "hg19/";
      const url = "https://genebe.net/variant/";
      let aa_regex = /^[ACGTacgt]+$/;
      if ( ref.match(aa_regex) && alt.match(aa_regex)) {
        let descriptor = `${chr}-${pos}-${ref}-${alt}`
        return `${url}${build}${descriptor}`;
      }
    }
    return ``
    }

    function genebe_viewer_link(row) {
      let ref = row["reference allele"];
      let alt = row["alternative allele"];
      if ( alt && ref) {
        let chr = row.chromosome;
        let pos = row.position;
        let gene = row.hgvsc.split('.')[0];
        const url = "https://genebe.net/tools/mutation-effect-viewer/";
        let aa_regex = /^[ACGTacgt]+$/;
        if ( ref.match(aa_regex) && alt.match(aa_regex)) {
          let descriptor = `${chr}-${pos}-${ref}-${alt}`
          return `${url}${gene}?mutations=${descriptor}`;
        }
      }
      return ``
    }

    function decipther_genomics_link(row) {
        gene = row.gene;
        return `https://www.deciphergenomics.org/gene/${gene}/overview/management-therapies/therapies`
    }

  function genomenexus_link(row) {
      let hgvsg = row.hgvsg
      if ( hgvsg ) {
        let build = ( row.build === "GRCh38" ) ? "grch38." : "";
        const url_suffix = "genomenexus.org/variant/"
        return `https://${build}${url_suffix}${hgvsg}`
      }
      return ``
    }

  function alphamissense_link(row) {
      if ( row.swissprot && row["protein position"] ) {
        let swissprot = row.swissprot.split('.')[0]
        let residue = row["protein position"]
        const url_prefix = "https://alphamissense.hegelab.org/hotspot?uid="
        return `${url_prefix}${swissprot}&resi=${residue}`
      }
      return ``
    }

  let hgvs_regex = /^(.+?):(.+)$/;
  // discard protein id
  let hgvsc = ( value ) ? value.match(hgvs_regex)[2].replace('%3D', '=') : "-";
  let links = "";


  let genebe_url = genebe_link(row);
  if (genebe_url != "") {
      links = `<a class="dropdown-item" href='${genebe_url}' target='_blank' rel='noopener noreferrer' >GeneBe</a>`
  }

  let genebe_viewer_url = genebe_viewer_link(row);
  if (genebe_viewer_url != "") {
      links = `${links}<a class="dropdown-item" href='${genebe_viewer_url}' target='_blank' rel='noopener noreferrer' >GeneBe Viewer</a>`
  }

  let decipther_genomics_url = decipther_genomics_link(row);
  links = `${links}<a class="dropdown-item" href='${decipther_genomics_url}' target='_blank' rel='noopener noreferrer' >Decipher Genomics</a>`

  let varsome_url = varsome_link(row);
  if (varsome_url != "") {
      links = `${links}<a class="dropdown-item" href='${varsome_url}' target='_blank' rel='noopener noreferrer' >VarSome</a>`
  }

  let genomenexus_url = genomenexus_link(row)
  if (genomenexus_url != "") {
      links = `${links}<a class="dropdown-item" href='${genomenexus_url}' target='_blank' rel='noopener noreferrer' >Genome Nexus</a>`
  }

  let alphamissense_url = alphamissense_link(row);
  if (alphamissense_url != "") {
      links = `${links}<a class="dropdown-item" href='${alphamissense_url}' target='_blank' rel='noopener noreferrer' >AlphaMissense</a>`
  }
  if ( links ) {
    let dropdown = `<div class="btn-group">
    <button class="btn btn-outline-secondary btn-table btn-sm dropdown-toggle" type="button" data-toggle="dropdown" aria-haspopup="true" aria-expanded="false">
    ${hgvsc}
    </button>
    <div class="dropdown-menu">
    ${links}
    </div>
    </div>
    `;
    return dropdown;
  }
  return ``
}
