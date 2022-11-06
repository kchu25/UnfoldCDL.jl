const html_template_has_alpha_valid=mt"""<!DOCTYPE html>
<html>
<meta charset="utf-8">
<head>
    <title></title>
    <style>
    abbr[title] {
    text-decoration: none;
    }
    table, td {
        border-collapse: collapse;
        margin: 15px 15px;
        padding: 5px 5px;
        table-layout: fixed;
        min-width: 85px;
    }
    .top_row {
    font-weight: bold;
    color: #808080;
    }

    thead,tfoot {
        font-weight: bold;
        background-color: #333;
        color:white;
    }

    .info {
        background-color: #E2E2E2;
        margin:5px;
        padding:5px;
    }
    </style>
    <script id="MathJax-script" async
        src="https://cdn.jsdelivr.net/npm/mathjax@3.0.0/es5/tex-mml-chtml.js">
</script>
</head>

<body>
    <div style="display:flex;">
        <div style="float:left; margin:25px; border:1px solid black; max-width:500px; padding:10px;" >
            Number of sequences: {{:num_seq}} <br>
            <table>
                <thead>
                    <tr>
                        <th colspan="100%">
                            Discovered motifs
                        </th>
                    </tr>
                </thead>
            <tbody>	
                    <tr class="top_row">
                        <td><abbr title="A label assigned to each discovered motif">Label</abbr></td>
                        <td><abbr title="P-value of the Fisher exact test on the odds ratio of the enrichment on the test set and the control dataset. The control dataset is the shuffled version of the test dataset that preserves the frequency of 2-mers in each seqeunce.">
                                        P-value</abbr></td>
                        <td><abbr title="Position weight matrix">Logo</abbr></td>
                        <td><abbr title="An estimate of the number of instances in the dataset"># instances</abbr></td>
                    </tr>		

                    {{#:DF}}
                    <tr>
                        <td style="text-align:center"><a href="{{:logo_folder}}/{{:logo}}.transfac">{{:label}}</a></td>
                        <td>{{{:eval}}}</td>
                        <td><img id="d_logo_{{:label}}" width="165" src="{{:logo_folder}}/{{:logo}}.png"><br>
                            <div id="d_orientation_{{:label}}">Learned PWM</div><br>
                            <button type="button" onclick="discovered_{{:label}}_changeToRC()">Reverse complement</button>
                        </td>
                        <td style="text-align:center">{{{:counts}}}</td>
                        <script type="text/javascript">					
                                    function discovered_{{:label}}_changeToRC() {
                                        var image = document.getElementById("d_logo_{{:label}}");
                                        if (image.src.match("_c")) {
                                            image.src = "{{:logo_folder}}/{{:logo}}.png";
                                        } else {
                                            image.src = "{{:logo_folder}}/{{:logo}}_c.png";
                                        }
                                        var orientation = document.getElementById("d_orientation_{{:label}}");
                                        if (orientation.innerHTML === "Learned PWM"){
                                            orientation.innerHTML = "Learned PWM's reverse-complement";
                                        } else {
                                            orientation.innerHTML = "Learned PWM";
                                        }
                                    }						 
                        </script>	
                    </tr>
                    {{/:DF}}		
                </tbody>
            </table>
            <br><br>	
        </div>	
        <div style="float:right; max-width:875px; margin:25px; padding:10px;"> 	
            <table>
                <thead>
                    <tr>
                        <th colspan="100%">
                            Motif co-occurrences
                        </th>
                    </tr>
                </thead>
                <tr><td style="text-align:center"> 
                    Context width: <input id="valR" type="range" min="0" max="{{:num_alphas}}" 
                        value="0" oninput="showVal(this.value)" onchange="showVal(this.value)"/>
                    <span id="range">{{:min_alpha}}</span>
                </td></tr>     
                <tr><td style="text-align:center"> <img id="img" width="100%" height="auto"></td></tr>                           
            </table>
            <script>
                var alphas ={{:valid_alphas}};
                var val = document.getElementById("valR").value;
                document.getElementById("range").innerHTML=alphas[val];
                document.getElementById("img").src = "logos/coocc_" + String(alphas[val]) + "_ww.png";
                function showVal(newVal){
                    document.getElementById("range").innerHTML=alphas[newVal];
                    document.getElementById("img").src = "logos/coocc_" + String(alphas[newVal]) + "_ww.png";
                }              
            </script>
        </div>
    </div>
</body>
</html>
"""

const html_template_no_alpha_valid=mt"""<!DOCTYPE html>
<html>
<meta charset="utf-8">
<head>
    <title></title>
    <style>
    abbr[title] {
    text-decoration: none;
    }
    table, td {
        border-collapse: collapse;
        margin: 15px 15px;
        padding: 5px 5px;
        table-layout: fixed;
        min-width: 85px;
    }
    .top_row {
    font-weight: bold;
    color: #808080;
    }

    thead,tfoot {
        font-weight: bold;
        background-color: #333;
        color:white;
    }

    .info {
        background-color: #E2E2E2;
        margin:5px;
        padding:5px;
    }
    </style>
    <script id="MathJax-script" async
        src="https://cdn.jsdelivr.net/npm/mathjax@3.0.0/es5/tex-mml-chtml.js">
</script>
</head>

<body>
    <div style="display:flex;">
        <div style="float:left; margin:25px; border:1px solid black; max-width:500px; padding:10px;" > 			
            Number of sequences: {{:num_seq}} <br>
            <table>
                <thead>
                    <tr>
                        <th colspan="100%">
                            Discovered motifs
                        </th>
                    </tr>
                </thead>
            <tbody>	
                    <tr class="top_row">
                        <td><abbr title="A label assigned to each discovered motif">Label</abbr></td>
                        <td><abbr title="P-value of the Fisher exact test on the odds ratio of the enrichment on the test set and the control dataset. The control dataset is the shuffled version of the test dataset that preserves the frequency of 2-mers in each seqeunce.">
                                        P-value</abbr></td>
                        <td><abbr title="Position weight matrix">Logo</abbr></td>
                        <td><abbr title="An estimate of the number of instances in the dataset"># instances</abbr></td>
                    </tr>		

                    {{#:DF}}
                    <tr>
                        <td style="text-align:center"><a href="{{:logo_folder}}/{{:logo}}.transfac">{{:label}}</a></td>
                        <td>{{{:eval}}}</td>
                        <td><img id="d_logo_{{:label}}" width="165" src="{{:logo_folder}}/{{:logo}}.png"><br>
                            <div id="d_orientation_{{:label}}">Learned PWM</div><br>
                            <button type="button" onclick="discovered_{{:label}}_changeToRC()">Reverse complement</button>
                        </td>
                        <td style="text-align:center">{{{:counts}}}</td>
                        <script type="text/javascript">					
                                    function discovered_{{:label}}_changeToRC() {
                                        var image = document.getElementById("d_logo_{{:label}}");
                                        if (image.src.match("_c")) {
                                            image.src = "{{:logo_folder}}/{{:logo}}.png";
                                        } else {
                                            image.src = "{{:logo_folder}}/{{:logo}}_c.png";
                                        }
                                        var orientation = document.getElementById("d_orientation_{{:label}}");
                                        if (orientation.innerHTML === "Learned PWM"){
                                            orientation.innerHTML = "Learned PWM's reverse-complement";
                                        } else {
                                            orientation.innerHTML = "Learned PWM";
                                        }
                                    }
                        </script>
                    </tr>
                    {{/:DF}}
                </tbody>
            </table>
            <br><br>
        </div>
    </div>
</body>
</html>
"""