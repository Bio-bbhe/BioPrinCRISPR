<template>
  <div id="main-container">
    <!-- Hero Section -->
    <header style="display: flex; flex-direction: row; min-height: 100vh;">
      <img :src="'/cover.png'" alt="Cas Protein Architecture Atlas" style="padding: 30px;padding-top: 0px;width:45%;height: 100vh;object-fit: contain;"/>
      <div class="hero">
        <h1>Cas Protein Architecture Atlas</h1>
        <h3>Mapping the Hidden Universe of CRISPR-Cas Systems</h3>
        <p>
            An interactive atlas of novel CRISPR-Cas systems discovered across 1,067,859 prokaryotic genomes. 
            Explore novel effector architectures, unique domain fusions, and the next generation of potential genome editing tools.
        </p>
        <a href="/#/network" target="_blank" class="btn btn-primary">Explore the Atlas →</a>
        <div class="secondary-links">
            <a href="/#/network/info?id=4830&label=SLFN_AlbA_2&example=1" target="_blank" class="btn btn-secondary">Read the Research Paper ↗</a>
            <a href="https://github.com/Bio-bbhe/BioPrinCRISPR" target="_blank" class="btn btn-secondary" >View Code on GitHub ↗</a>
        </div>
      </div>
    </header>

    <main>
      <!-- Principles Section -->
      <section id="principles" class="container">
        <h2>A Discovery Engine Built on First Principles</h2>
        <p class="section-subtitle">
          Our platform integrates multiple layers of biological data to move beyond conventional homology search. 
          We connect co-conservation, domain architecture, and deep learning embeddings to reveal systems that were previously invisible.
        </p>
        <div class="principles-grid">
          <div class="principle-item">
            <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><path d="M12 2C6.486 2 2 6.486 2 12s4.486 10 10 10 10-4.486 10-10S17.514 2 12 2zm0 18c-4.411 0-8-3.589-8-8s3.589-8 8-8 8 3.589 8 8-3.589 8-8 8z"></path><path d="M11 11h2v2h-2zm-4 0h2v2H7zm8 0h2v2h-2z"></path></svg>
            <h4>Domain Networks</h4>
            <p>Visualize the universe of protein domain relationships and uncover novel functional linkages.</p>
          </div>
          <div class="principle-item">
            <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><path d="M11.293 2.293 2.586 11H8v9h8v-9h5.414L12.707 2.293a1 1 0 0 0-1.414 0zM18 20v-8h1.414L12 4.586 4.586 12H6v8h12z"></path></svg>
            <h4>Genomic Context</h4>
            <p>Browse interactive operon maps for every system to understand gene organization and co-conservation.</p>
          </div>
          <div class="principle-item">
            <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24"><path d="m21.446 11.23-8-5.5a1 1 0 0 0-1.012.042l-8 6A1 1 0 0 0 4 12.618V17a1 1 0 0 0 1 1h14a1 1 0 0 0 1-1v-4.5a1 1 0 0 0-.554-.882zM18 16H6v-3.382l6-4.5 6 4.364V16z"></path><path d="M12 18.5a1 1 0 0 0-1 1V21a1 1 0 0 0 2 0v-1.5a1 1 0 0 0-1-1z"></path></svg>
            <h4>3D Protein Structures</h4>
            <p>Inspect predicted 3D structures for novel effectors to gain insights into their potential mechanisms.</p>
          </div>
        </div>
      </section>

      <!-- Discoveries Section -->
      <section id="discoveries" class="container discoveries">
        <h2 style="color: #ffff;">Spotlight on Novel Discoveries</h2>
        <p class="section-subtitle">
          A curated look at some of the most significant findings from the BioPrinCRISPR pipeline, now available for exploration in the Atlas.
        </p>
        <div class="discoveries-grid">
          <div class="discovery-card" v-for="(item, i) in discoveries" :key="i">
            <!-- <img :src="item.img" :alt="item.alt" /> -->
            <div class="card-content">
              <h4>{{ item.title }}</h4>
              <p>{{ item.desc }}</p>
              <a :href="item.link" target="_blank">{{ item.linkText }}</a>
            </div>
          </div>
        </div>
      </section>
    </main>

    <!-- Footer -->
    <footer>
      <div class="container footer-content">
        <svg viewBox="0 0 24 24"><path d="M12 2.5a3.5 3.5 0 0 0-3.5 3.5v.5h7v-.5a3.5 3.5 0 0 0-3.5-3.5zM8.5 7v1.126a3.5 3.5 0 0 0-1.5 2.874v5a1 1 0 0 0 1 1h7a1 1 0 0 0 1-1v-5a3.5 3.5 0 0 0-1.5-2.874V7h-5z"></path></svg>
        <p>Developed in the [Your Lab Name] at [Your University]</p>
        <p>Funding by [Your Funding Agency]</p>
      </div>
    </footer>
  </div>
</template>

<script setup>
import { onMounted } from 'vue';

const discoveries = [
    {
        img: "https://via.placeholder.com/400x180/000000/7E57C2.png?text=EpaF+Structure",
        alt: "EpaF Structure",
        title: "A New Class of Effector: The EpaF Nuclease-Helicase",
        desc: "Discover a unique single-protein effector that fuses DNA unwinding and cleavage domains, representing a novel architecture for CRISPR-mediated defense.",
        link: "/#/network/info?id=3927&label=PDDEXK_1&example=1",
        linkText: "Explore EpaF Systems →"
    },
    {
        img: "https://via.placeholder.com/400x180/000000/7E57C2.png?text=Validated+RNA+Editing",
        alt: "Cas3–Cse1 Fusions",
        title: "Non-canonical Effectors: Cas3–Cse1 Fusions",
        desc: "Highlighting atypical hybrid proteins where the Cas3 C-terminal region fuses with Cse1, forming non-canonical Class 1 effectors that may enable specialized CRISPR interference strategies.",
        link: "/#/network/info?id=888&label=Cas3_HD&example=1",
        linkText: "Explore Cas3 Systems →"
    },
    {
        img: "https://via.placeholder.com/400x180/000000/7E57C2.png?text=Accessory+Module",
        alt: "Accessory Module",
        title: "Novel Accessory Modules: Schlafen-Res Fusions",
        desc: "Investigate a chimeric protein combining helicase and nuclease domains that may enhance canonical interference, found in Type III and Class I loci.",
        link: "/#/network/info?id=6039&label=SLFN-g3_helicase&example=1",
        linkText: "Explore Accessory Modules →"
    }
];

onMounted(() => {
    document.title = 'Cas Protein Architecture Atlas';
});
</script>

<style scoped>
/* --- General Setup & Fonts --- */
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;700&display=swap');

* {
    margin: 0;
    padding: 0;
    box-sizing: border-box;
}

#main-container {
    font-family: 'Inter', sans-serif;
    background-color: var(--color-background);
    color: var(--color-text);
    line-height: 1.6;
    overflow-x: hidden;

    --color-background: #ffffff;
    --color-text: #000000;
    --color-text-muted: #888;
    --color-primary: #ffcf25;
    --color-primary-glow: rgba(175, 175, 175, 0.5);
    --color-surface: #1a1a2e;
}

.particle {
    position: absolute;
    border-radius: 50%;
    background: var(--color-text);
    opacity: 0;
    animation: float 20s infinite linear;
}

@keyframes float {
    0% { transform: translateY(0) translateX(0); opacity: 0; }
    10% { opacity: 0.8; }
    90% { opacity: 0.8; }
    100% { transform: translateY(-100vh) translateX(5vw); opacity: 0; }
}

/* --- Reusable Components --- */
.container {
    max-width: 1100px;
    margin: 0 auto;
    padding: 4rem 2rem;
}

.btn {
    display: inline-block;
    padding: 12px 28px;
    border-radius: 8px;
    text-decoration: none;
    font-weight: 500;
    transition: all 0.3s ease;
}

.btn-primary {
    background-color: var(--color-primary);
    color: white;
    font-size: 1.1rem;
    box-shadow: 0 4px 20px var(--color-primary-glow);
}

.btn-primary:hover {
    transform: translateY(-3px);
    box-shadow: 0 6px 25px var(--color-primary-glow);
}

.btn-secondary {
    border: 1px solid var(--color-text-muted);
    color: var(--color-text);
    margin: 0 0.5rem;
}

.btn-secondary:hover {
    background-color: var(--color-surface);
    border-color: var(--color-text);
    color: white;
}

/* --- Header / Hero Section --- */
.hero {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    text-align: center;
    padding: 0 1rem;
    width:55%;
}

.hero h1 {
    font-size: 3.5rem;
    font-weight: 700;
    margin-bottom: 0.5rem;
    letter-spacing: -1px;
}

.hero h3 {
    font-size: 1.5rem;
    font-weight: 400;
    margin-bottom: 1.5rem;
    color: var(--color-text-muted);
}

.hero p {
    max-width: 600px;
    margin-bottom: 2.5rem;
}

.hero .secondary-links {
    margin-top: 2rem;
}

/* --- Sections Styling --- */
section {
    text-align: center;
}

section h2 {
    font-size: 2.5rem;
    margin-bottom: 1rem;
}

section .section-subtitle {
    max-width: 700px;
    margin: 0 auto 3rem auto;
    color: var(--color-text-muted);
}

/* --- Principles Section --- */
.principles-grid {
    display: flex;
    justify-content: space-around;
    gap: 2rem;
    margin-top: 3rem;
}

.principle-item {
    flex: 1;
    padding: 1.5rem;
}

.principle-item svg {
    width: 50px;
    height: 50px;
    margin-bottom: 1rem;
    fill: var(--color-primary);
}

.principle-item h4 {
    font-size: 1.2rem;
    margin-bottom: 0.5rem;
}

.principle-item p {
    color: var(--color-text-muted);
    font-size: 0.9rem;
}

/* --- Discoveries Section --- */
.discoveries {
    background-color: var(--color-surface);
    border-radius: 20px;
}

.discoveries-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
    gap: 2rem;
}

.discovery-card {
    background-color: var(--color-background);
    border-radius: 12px;
    overflow: hidden;
    text-align: left;
    transition: transform 0.3s ease, box-shadow 0.3s ease;
}

.discovery-card:hover {
    transform: translateY(-8px);
    box-shadow: 0 10px 30px rgba(0, 0, 0, 0.5);
}

.discovery-card img {
    width: 100%;
    height: 180px;
    object-fit: cover;
    background-color: #000;
}

.card-content {
    padding: 1.5rem;
}

.card-content h4 {
    font-size: 1.3rem;
    margin-bottom: 0.5rem;
}

.card-content p {
    font-size: 0.95rem;
    color: var(--color-text-muted);
    margin-bottom: 1.5rem;
}

.card-content a {
    color: var(--color-primary);
    text-decoration: none;
    font-weight: 500;
}

/* --- Footer --- */
footer {
    padding: 3rem 2rem;
    text-align: center;
}

.footer-content {
    display: flex;
    justify-content: space-between;
    align-items: center;
    max-width: 800px;
    margin: 0 auto;
    opacity: 0.7;
    font-size: 0.9rem;
}

.footer-content svg {
    height: 40px;
    fill: var(--color-text-muted);
}

/* --- Responsive Design --- */
@media (max-width: 768px) {
    .hero h1 { font-size: 2.5rem; }
    .hero h3 { font-size: 1.2rem; }
    .principles-grid { flex-direction: column; }
    .footer-content { flex-direction: column; gap: 1rem; }
}
</style>
